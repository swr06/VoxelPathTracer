#version 450 core

#define CLOUD_HEIGHT 70
#define PCF_COUNT 4 // we really dont care about the low sample count because we have taa.
//#define POISSON_DISK_SAMPLING
#define PI 3.14159265359
#define THRESH 1.41414


layout (location = 0) out vec3 o_Color;
layout (location = 1) out vec3 o_Normal;
layout (location = 2) out vec4 o_PBR;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_DiffuseTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InitialTracePositionTexture;
uniform sampler2D u_BlockIDTexture;
uniform sampler2D u_ShadowTexture;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlueNoiseTextures;
uniform sampler2DArray u_BlockEmissiveTextures;
uniform samplerCube u_Skybox;
uniform sampler2D u_CloudData;
uniform sampler2D u_PreviousNormalTexture; 
uniform sampler2D u_VXAO;

uniform bool u_UseDFG;
uniform bool u_RemoveTiling;

uniform sampler2D u_DiffuseSHData1;
uniform sampler2D u_DiffuseSHData2;

uniform sampler2D u_ReflectionSHData;
uniform sampler2D u_ReflectionCoCgData;

uniform sampler2D u_HighResBL;

uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;
uniform vec3 u_StrongerLightDirection;

uniform float u_Time;
uniform float u_GrassblockAlbedoID;

uniform mat4 u_ShadowView;
uniform mat4 u_ShadowProjection;
uniform mat4 u_ReflectionView;
uniform mat4 u_ReflectionProjection;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform vec3 u_ViewerPosition;
uniform vec2 u_Dimensions;

uniform bool u_CloudsEnabled;
uniform bool u_POM = false;
uniform bool u_HighQualityPOM = false;
uniform bool u_RTAO;
uniform bool u_ContactHardeningShadows;
uniform bool u_AmplifyNormalMap;

uniform bool u_DoVXAO = true;
uniform bool u_SVGFEnabled = true;
uniform bool u_ShouldDitherUpscale = true;

uniform float u_CloudBoxSize;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

// *****Temporary****** solution to have multi texturing for grass blocks
// Data stored : 
// Grass block ID
// top face index (albedo, normal, pbr),
// right/left/front/back face index (albedo, normal, pbr),
// bottom face index (albedo, normal, pbr)
// 9 + 1 ints total

uniform int u_GrassBlockProps[10];


vec4 textureBicubic(sampler2D sampler, vec2 texCoords);
vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 radiance_s, vec3 albedo, vec3 normal, vec3 pbr, float shadow);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);


const vec3 ATMOSPHERE_SUN_COLOR = vec3(1.0f, 1.0f, 0.5f);
const vec3 ATMOSPHERE_MOON_COLOR =  vec3(0.1f, 0.1f, 1.0f);


vec3 PoissonDisk3D[16] = vec3[](
    vec3(0.488937, 0.374798, 0.0314035),
    vec3(0.112522, 0.0911893, 0.932066),
    vec3(0.345347, 0.857173, 0.622028),
    vec3(0.724845, 0.0422376, 0.754479),
    vec3(0.775262, 0.82693, 0.16596),
    vec3(0.971221, 0.020539, 0.0113529),
    vec3(0.010834, 0.171209, 0.379254),
    vec3(0.577593, 0.514908, 0.977874),
    vec3(0.170507, 0.840266, 0.0510269),
    vec3(0.939055, 0.566179, 0.568987),
    vec3(0.015137, 0.606647, 0.998566),
    vec3(0.687002, 0.465712, 0.479293),
    vec3(0.170232, 0.25837, 0.602069),
    vec3(0.83755 , 0.334819, 0.0497452),
    vec3(0.795679, 0.742149, 0.878201),
    vec3(0.180761, 0.585253, 0.245888)
);


float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}
#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))

const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

bool CompareFloatNormal(float x, float y) {
    return abs(x - y) < 0.02f;
}

vec3 GetNormalFromID(float n) {
	const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f));
    int idx = int(round(n*10.0f));

    if (idx > 5) {
        return vec3(1.0f, 1.0f, 1.0f);
    }

    return Normals[idx];
}

vec3 SampleNormal(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}


float Noise2d( in vec2 x )
{
    float xhash = cos( x.x * 37.0 );
    float yhash = cos( x.y * 57.0 );
    return fract( 415.92653 * ( xhash + yhash ) );
}

// thresholded white noise :
float NoisyStarField( in vec2 vSamplePos, float fThreshhold )
{
    float StarVal = Noise2d(vSamplePos);
    if ( StarVal >= fThreshhold )
        StarVal = pow((StarVal - fThreshhold)/(1.0 - fThreshhold), 6.0 );
    else
        StarVal = 0.0;
    return StarVal;
}

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SamplePositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

float StableStarField(in vec2 vSamplePos, float fThreshhold)
{
    float fractX = fract(vSamplePos.x);
    float fractY = fract(vSamplePos.y);
    vec2 floorSample = floor( vSamplePos );
    float v1 = NoisyStarField( floorSample, fThreshhold );
    float v2 = NoisyStarField( floorSample + vec2(0.0, 1.0), fThreshhold );
    float v3 = NoisyStarField( floorSample + vec2(1.0, 0.0), fThreshhold );
    float v4 = NoisyStarField( floorSample + vec2(1.0, 1.0), fThreshhold );
    float StarVal =   v1 * (1.0 - fractX) * (1.0 - fractY)
        			+ v2 * (1.0 - fractX) * fractY
        			+ v3 * fractX * (1.0 - fractY)
        			+ v4 * fractX * fractY;
	return StarVal;
}

float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float stars(vec3 fragpos)
{
    fragpos.y = abs(fragpos.y);
	float elevation = clamp(fragpos.y, 0.0f, 1.0f);
	vec2 uv = fragpos.xz / (1.0f + elevation);
    float star = StableStarField(uv * 700.0f, 0.98);
    star *= 0.06f;
    float rand_val = rand(fragpos.xy);
	return clamp(star, 0.0f, 100000.0f) * 4.46f;
}

vec2 RayBoxIntersect(vec3 boundsMin, vec3 boundsMax, vec3 rayOrigin, vec3 invRaydir)
{
	vec3 t0 = (boundsMin - rayOrigin) * invRaydir;
	vec3 t1 = (boundsMax - rayOrigin) * invRaydir;
	vec3 tmin = min(t0, t1);
	vec3 tmax = max(t0, t1);
	float dstA = max(max(tmin.x, tmin.y), tmin.z);
	float dstB = min(tmax.x, min(tmax.y, tmax.z));
	float dstToBox = max(0, dstA);
	float dstInsideBox = max(0, dstB - dstToBox);
	return vec2(dstToBox, dstInsideBox);
}

bool GetAtmosphere(inout vec3 atmosphere_color, in vec3 in_ray_dir, float transmittance, float true_transmittance)
{
    vec3 sun_dir = normalize(u_SunDirection); 
    vec3 moon_dir = vec3(-sun_dir.x, -sun_dir.y, sun_dir.z); 

    vec3 ray_dir = normalize(in_ray_dir);
    
    if (transmittance > 0.7f) {

        if(dot(ray_dir, sun_dir) > 0.999825f)
        {
            atmosphere_color = ATMOSPHERE_SUN_COLOR; 
            o_PBR.w = float(2.6f);
            return true;
        }
        
        if(dot(ray_dir, moon_dir) > 0.99986f)
        {
            atmosphere_color = ATMOSPHERE_MOON_COLOR;
            o_PBR.w = float(1.2f);
            return true;
        }
    }

    vec3 atmosphere = texture(u_Skybox, ray_dir).rgb;

    float star_visibility;
    star_visibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; 
    star_visibility = 1.0f - star_visibility;
    vec3 stars = vec3(stars(vec3(in_ray_dir)) * star_visibility);
    stars = clamp(stars, 0.0f, 1.3f);

    atmosphere += stars*4.0f*true_transmittance;

    atmosphere_color = atmosphere;

    return false;
}

float CASWeight(vec3 x) {
    //return GetSat(x);
    return max(dot(x, vec3(0.2125f, 0.7154f, 0.0721f)), 0.01f);
    //return max(x.g, EPSILON);
}

vec4 CAS(sampler2D Texture, vec2 uv, float SharpeningAmount)
{
    ivec2 Pixel = ivec2(floor(uv*textureSize(Texture,0)));

    // Samples 
    vec4 a = texelFetch(Texture, Pixel + ivec2(0, -1), 0).rgba;
    vec4 b = texelFetch(Texture, Pixel + ivec2(-1, 0), 0).rgba;
    vec4 c = texelFetch(Texture, Pixel + ivec2(0, 0), 0).rgba;
    vec4 d = texelFetch(Texture, Pixel + ivec2(1, 0), 0).rgba;
    vec4 e = texelFetch(Texture, Pixel + ivec2(0, 1), 0).rgba;

    // Weight by luminance 
    float WeightA = CASWeight(a.xyz);
    float WeightB = CASWeight(b.xyz);
    float WeightC = CASWeight(c.xyz);
    float WeightD = CASWeight(d.xyz);
    float WeightE = CASWeight(e.xyz);

    // Calculate bounds :
    float MinWeighter = min(WeightA, min(WeightB, min(WeightC, min(WeightD, WeightE))));
    float MaxWeighter = max(WeightA, max(WeightB, max(WeightC, max(WeightD, WeightE))));

    // Apply weights :
    float FinalSharpenAmount = sqrt(min(1.0f - MaxWeighter, MinWeighter) / MaxWeighter);
    float w = FinalSharpenAmount * mix(-0.125f, -0.2f, SharpeningAmount);
    return (w * (a + b + d + e) + c) / (4.0f * w + 1.0f);
}

vec3 GetAtmosphereAndClouds()
{
    vec2 ij = floor(mod(gl_FragCoord.xy, vec2(2.0) ));
	float idx = ij.x + 2.0*ij.y;
	vec4 m = step( abs(vec4(idx)-vec4(0,1,2,3)), vec4(0.5) ) * vec4(0.75,0.25,0.00,0.50);
	float d = m.x+m.y+m.z+m.w; // dither.
    float base_bayer = bayer4(gl_FragCoord.xy);
    vec3 NormalizedDir = normalize(v_RayDirection);
    float CloudFade  = mix(0.0f, 1.0f, max(NormalizedDir.y, 0.00001f));
    CloudFade = pow(CloudFade * 2.75f, 3.25f);
    CloudFade = clamp(CloudFade, 0.0f, 1.0f);
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
    const vec3 D = (vec3(355.0f, 10.0f, 0.0f) / 255.0f) * 0.4f;
    vec3 S = vec3(1.45f);
    float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 1.8f), 0.0f, 1.0f);
    S = mix(S, D, DuskVisibility);
    vec3 M = mix(S + 0.001f, (vec3(46.0f, 142.0f, 255.0f) / 255.0f) * 0.2f, SunVisibility); 
	vec4 SampledCloudData = texture(u_CloudData, v_TexCoords+(base_bayer*(1.0f/textureSize(u_CloudData,0)))).rgba; // Bicubic B spline interp
	//vec4 SampledCloudData = texture(u_CloudData, v_TexCoords).rgba; // Bicubic B spline interp
	//vec4 SampledCloudData = CAS(u_CloudData,v_TexCoords,0.8f);
    SampledCloudData.xyz *= CloudFade;
    float Transmittance = SampledCloudData.w;
    vec3 Scatter = SampledCloudData.xyz*1.2f ;
    vec3 Sky = vec3(0.0f);
    bool v = GetAtmosphere(Sky, NormalizedDir, Transmittance*6.5f, Transmittance);
    Sky = Sky * max(Transmittance, 0.9f);
    return vec3(Sky + Scatter*M)+(d/255.0f); // + dither.
}

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

vec4 ClampedTexture(sampler2D tex, vec2 txc)
{
    return texture(tex, clamp(txc, 0.0f, 1.0f));
}

vec2 ReprojectShadow(in vec3 world_pos)
{
	vec3 WorldPos = world_pos;

	vec4 ProjectedPosition = u_ShadowProjection * u_ShadowView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

int MIN = -2147483648;
int MAX = 2147483647;

int xorshift(in int value) 
{
    // Xorshift*32
    // Based on George Marsaglia's work: http://www.jstatsoft.org/v08/i14/paper
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    return value;
}

int nextInt(inout int seed) 
{
    seed = xorshift(seed);
    return seed;
}

float nextFloat(inout int seed) 
{
    seed = xorshift(seed);
    // FIXME: This should have been a seed mapped from MIN..MAX to 0..1 instead
    return abs(fract(float(seed) / 3141.592653));
}

float nextFloat(inout int seed, in float max) 
{
    return nextFloat(seed) * max;
}

float nextFloat(inout int seed, in float min, in float max) 
{
    return min + (max - min) * nextFloat(seed);
}


bool RayBoxIntersect(const vec3 boxMin, const vec3 boxMax, vec3 r0, vec3 rD, out float t_min, out float t_max) 
{
	vec3 inv_dir = 1.0f / rD;
	vec3 tbot = inv_dir * (boxMin - r0);
	vec3 ttop = inv_dir * (boxMax - r0);
	vec3 tmin = min(ttop, tbot);
	vec3 tmax = max(ttop, tbot);
	vec2 t = max(tmin.xx, tmin.yz);
	float t0 = max(t.x, t.y);
	t = min(tmax.xx, tmax.yz);
	float t1 = min(t.x, t.y);
	t_min = t0;
	t_max = t1;
	return t1 > max(t0, 0.0);
}

// http://www.iquilezles.org/www/articles/intersectors/intersectors.htm
float capIntersect( in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float r )
{
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;

    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoa = dot(ba,oa);
    float rdoa = dot(rd,oa);
    float oaoa = dot(oa,oa);

    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - r*r*baba;
    float h = b*b - a*c;
    if( h>=0.0 )
    {
        float t = (-b-sqrt(h))/a;
        float y = baoa + t*bard;
        if( y>0.0 && y<baba ) return t;
        vec3 oc = (y<=0.0) ? oa : ro - pb;
        b = dot(rd,oc);
        c = dot(oc,oc) - r*r;
        h = b*b - c;
        if( h>0.0 ) return -b - sqrt(h);
    }
    return -1.0;
}

bool GetPlayerIntersect(in vec3 WorldPos, in vec3 d)
{
    float t = capIntersect(WorldPos, d, u_ViewerPosition, u_ViewerPosition + vec3(0.0f, 1.0f, 0.0f), 0.5f);
    return t > 0.0f;
}

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

bool InScreenSpace(vec2 x)
{
    return x.x < 1.0f && x.x > 0.0f && x.y < 1.0f && x.y > 0.0f;
}

int RNG_SEED;

float ComputeShadow(vec3 world_pos, vec3 flat_normal)
{
    vec2 Txc = ReprojectShadow(world_pos);
    float BaseShadow = texture(u_ShadowTexture, Txc).r;
    float Shadow = BaseShadow;
    float PlayerShadow  = 0.0f;
    int PlayerShadowSamples = 16;

    vec3 BiasedWorldPos = world_pos - (flat_normal * 0.255f);
    if (u_ContactHardeningShadows) {
        for (int i = 0 ; i < PlayerShadowSamples ; i++)
        {
            vec2 Hash = hash2();
	        vec2 Hash2 = hash2();
            vec3 JitteredLightDirection = u_StrongerLightDirection;
            JitteredLightDirection.x += Hash.x * 0.03;
	        JitteredLightDirection.z += Hash.y * 0.03;
	        JitteredLightDirection.y += abs(Hash2.y * 0.01);
	        JitteredLightDirection = normalize(JitteredLightDirection);
            PlayerShadow += float(GetPlayerIntersect(BiasedWorldPos.xyz, JitteredLightDirection.xyz));
        }

        PlayerShadow /= float(PlayerShadowSamples);
    }

    else 
    {
        PlayerShadow = float(GetPlayerIntersect(world_pos.xyz, normalize(u_StrongerLightDirection).xyz));
    }

    return clamp(BaseShadow + PlayerShadow, 0.0f, 1.0f);
}

bool IsInScreenSpaceBounds(in vec2 tx)
{
    if (tx.x > 0.0f && tx.y > 0.0f && tx.x < 1.0f && tx.y < 1.0f)
    {
        return true;
    }

    return false;
}

float GetDisplacementAt(in vec2 txc, in float pbridx) 
{
    return texture(u_BlockPBRTextures, vec3(vec2(txc.x, txc.y), pbridx)).b * 0.35f;
}

vec2 ParallaxOcclusionMapping(vec2 TextureCoords, vec3 ViewDirection, in float pbridx) // View direction should be in tangent space!
{ 
    float NumLayers = u_HighQualityPOM ? 72 : 32; 
    float LayerDepth = 1.0 / (u_HighQualityPOM ? (NumLayers * 0.725f) : (NumLayers * 0.65));
    float CurrentLayerDepth = 0.0;
    vec2 P = ViewDirection.xy * 1.0f; 
    vec2 DeltaTexCoords = P / NumLayers;
    vec2 InitialDeltaCoords = DeltaTexCoords;

    vec2  CurrentTexCoords = TextureCoords;
    float CurrentDepthMapValue = GetDisplacementAt(CurrentTexCoords, pbridx);

    for (int i = 0 ; i < NumLayers ; i++)
    {
        if(CurrentLayerDepth < CurrentDepthMapValue)
        {
            //CurrentTexCoords -= max(DeltaTexCoords, InitialDeltaCoords * 0.025f);
            CurrentTexCoords -= DeltaTexCoords;
            CurrentDepthMapValue = GetDisplacementAt(CurrentTexCoords, pbridx);  
            CurrentLayerDepth += LayerDepth;
            DeltaTexCoords *= u_HighQualityPOM ? 0.970f : 0.9f;
        }
    }

    vec2 PrevTexCoords = CurrentTexCoords + DeltaTexCoords;
    float AfterDepth  = CurrentDepthMapValue - CurrentLayerDepth;
    float BeforeDepth = GetDisplacementAt(PrevTexCoords, pbridx) - CurrentLayerDepth + LayerDepth;
    float Weight = AfterDepth / (AfterDepth - BeforeDepth);
    vec2 FinalTexCoords = PrevTexCoords * Weight + CurrentTexCoords * (1.0 - Weight);
    return FinalTexCoords;
}   

bool CompareVec3(vec3 v1, vec3 v2) {
	float e = 0.0125f;
	return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
}

void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );

	const vec3 Tangents[6] = vec3[]( vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(0.0f, 0.0f, -1.0f), vec3(0.0f, 0.0f, -1.0f)
				   );

	const vec3 BiTangents[6] = vec3[]( vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f),
				     vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 1.0f),
					 vec3(0.0f, -1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f)
	);

    uv = vec2(1.0f);

	if (CompareVec3(normal, Normals[0]))
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[0];
		bitangent = BiTangents[0];
    }

    else if (CompareVec3(normal, Normals[1]))
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[1];
		bitangent = BiTangents[1];
    }

    else if (CompareVec3(normal, Normals[2]))
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[2];
		bitangent = BiTangents[2];
    }

    else if (CompareVec3(normal, Normals[3]))
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[3];
		bitangent = BiTangents[3];
    }
	
    else if (CompareVec3(normal, Normals[4]))
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[4];
		bitangent = BiTangents[4];
    }
    

    else if (CompareVec3(normal, Normals[5]))
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[5];
		bitangent = BiTangents[5];
    }
}

vec3 fresnelroughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
{
    float cosTheta = max(dot(norm, Eye), 0.0);
    const float magic = 2.4f;
    return F0 + (max(vec3(pow(1.0f - roughness, magic)), F0) - F0) * pow(clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
}

vec4 GetTextureIDs(int BlockID) 
{
	return vec4(float(BlockAlbedoData[BlockID]),
				float(BlockNormalData[BlockID]),
				float(BlockPBRData[BlockID]),
				float(BlockEmissiveData[BlockID]));
}

int GetBlockID(vec2 txc)
{
	float id = texture(u_BlockIDTexture, txc).r;

	return clamp(int(floor(id * 255.0f)), 0, 127);
}

float FastDistance(vec3 p1, vec3 p2)
{
    return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

bool IsAtEdge(in vec2 txc)
{
    vec2 TexelSize = 1.0f / textureSize(u_InitialTracePositionTexture, 0);

    const vec2 Kernel[8] = vec2[8]
	(
		vec2(-1.0, -1.0),
		vec2( 0.0, -1.0),
		vec2( 1.0, -1.0),
		vec2(-1.0,  0.0),
		vec2( 1.0,  0.0),
		vec2(-1.0,  1.0),
		vec2( 0.0,  1.0),
		vec2( 1.0,  1.0)
	);

    for (int i = 0 ; i < 8 ; i++)
    {
        if (SamplePositionAt(u_InitialTracePositionTexture, txc + Kernel[i] * TexelSize).w <= 0.0f)
        {
            return true;
        }
    }

    return false;
}

#define SPATIAL_UPSCALE_INDIRECT_DIFFUSE

void SpatialUpscaleIndirectDiffuse(vec3 BaseNormal, float BaseLinearDepth, out vec4 SH, out vec2 CoCg)
{
#ifndef SPATIAL_UPSCALE_INDIRECT_DIFFUSE
    SH = texture(u_DiffuseSHData1, v_TexCoords).xyzw ;
	CoCg = texture(u_DiffuseSHData2, v_TexCoords).xy ;
    return;

#else 

	vec4 TotalSH = vec4(0.0f);
	vec2 TotalCoCg = vec2(0.0f);
	float TotalWeight = 0.0f;

	vec2 TexelSize = 1.0f / textureSize(u_DiffuseSHData1, 0);
	const float Weights[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);
    float HighFreqBL = texture(u_HighResBL, v_TexCoords * (u_Dimensions / vec2(1024.0f))).r;
    float SpatialDither = u_ShouldDitherUpscale ? bayer8(gl_FragCoord.xy) : 1.0f;

	for (int x = -2 ; x <= 2 ; x++) {
		for (int y = -2 ; y <= 2 ; y++) {
			
            vec2 SampleCoord = v_TexCoords + (vec2(x,y) * SpatialDither * 1.025f) * TexelSize;
			float LinearDepthAt = (texture(u_InitialTracePositionTexture, SampleCoord).x);
            if (LinearDepthAt <= 0.0f + 1e-2) { continue; }
			
            vec3 NormalAt = SampleNormal(u_NormalTexture, SampleCoord.xy).xyz;
			float DepthWeight = 1.0f / (abs(LinearDepthAt - BaseLinearDepth) + 0.01f);
			DepthWeight = pow(DepthWeight, 16.0f);
            DepthWeight = max(DepthWeight, 0.0f);
			float NormalWeight = pow(abs(dot(NormalAt, BaseNormal)), 8.0f);
            NormalWeight = max(NormalWeight, 0.0f);
			float Kernel = Weights[x + 2] * Weights[y + 2];
			float Weight = Kernel * clamp(NormalWeight * DepthWeight, 0.0f, 1.0f);
			Weight = max(Weight, 0.01f);
			TotalSH += texture(u_DiffuseSHData1, SampleCoord).xyzw * Weight;
			TotalCoCg += texture(u_DiffuseSHData2, SampleCoord).xy * Weight;
			TotalWeight += Weight;
		}
	}

    SH = TotalSH / TotalWeight;
    CoCg = TotalCoCg / TotalWeight;
#endif
}

vec3 saturate(vec3 x)
{
    return clamp(x, 0.0f, 1.0f);
}

// from quake2rtx
vec3 SHToIrridiance(vec4 shY, vec2 CoCg, vec3 v)
{
    float x = dot(shY.xyz, v);
    float Y = 2.0 * (1.023326f * x + 0.886226f * shY.w);
    Y = max(Y, 0.0);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

// from unreal
vec3 DFGPolynomialApproximate(vec3 F0, float NdotV, float roughness) 
{
    const vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
    const vec4 c1 = vec4(1.0, 0.0425, 1.04, -0.04);
    vec4 r = roughness * c0 + c1;
    float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
    vec2 AB = vec2(-1.04, 1.04) * a004 + r.zw;
    return F0 * AB.x + AB.y;
}

// COLORS //
const vec3 SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * 6.4f;
const vec3 NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.5f; 
const vec3 DUSK_COLOR = (vec3(255.0f, 204.0f, 144.0f) / 255.0f) * 0.064f; 


float sum( vec3 v ) { return v.x+v.y+v.z; }


vec4 hash4( vec2 p ) { return fract(sin(vec4( 1.0+dot(p,vec2(37.0,17.0)), 
                                              2.0+dot(p,vec2(11.0,47.0)),
                                              3.0+dot(p,vec2(41.0,29.0)),
                                              4.0+dot(p,vec2(23.0,31.0))))*103.0); }




void main()
{
	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(100.0f * fract(u_Time));

    // Xorshift!
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
    HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(u_Time) * 100.0f;

    vec4 WorldPosition = SamplePositionAt(u_InitialTracePositionTexture, v_TexCoords);

    vec3 SampledNormals = SampleNormal(u_NormalTexture, v_TexCoords).rgb;
    vec3 AtmosphereAt = vec3(0.0f);
    o_Color = vec3(1.0f);

    if (WorldPosition.w > 0.0f)
    {
        float RayTracedShadow = 0.0f;
        
        RayTracedShadow = ComputeShadow(WorldPosition.xyz, SampledNormals.xyz);

        float Bias = 0.0035f;
        bool InBiasedSS =  (v_TexCoords.x > 0.0 + Bias && v_TexCoords.x < 1.0 - Bias 
		 && v_TexCoords.y > 0.0 + Bias && v_TexCoords.y < 1.0f - Bias);

        if (!IsAtEdge(v_TexCoords))
        {
            vec2 UV;
            vec3 Tangent, Bitangent;

            CalculateVectors(WorldPosition.xyz, SampledNormals, Tangent, Bitangent, UV); 

	        mat3 tbn = mat3(normalize(Tangent), normalize(Bitangent), normalize(SampledNormals));
            int BaseBlockID = GetBlockID(v_TexCoords);
            vec4 data = GetTextureIDs(BaseBlockID);


            // Handle grass block! 
            // dont kill me pls :( 
            if (BaseBlockID == u_GrassBlockProps[0])
	        {
	            if (SampledNormals == NORMAL_LEFT || SampledNormals == NORMAL_RIGHT || SampledNormals == NORMAL_FRONT || SampledNormals == NORMAL_BACK)
	        	{
	        		data.x = u_GrassBlockProps[4];
	        		data.y = u_GrassBlockProps[5];
	        		data.z = u_GrassBlockProps[6];
	        	}

	        	else if (SampledNormals == NORMAL_TOP)
	        	{
	        		data.x = u_GrassBlockProps[1];
	        		data.y = u_GrassBlockProps[2];
	        		data.z = u_GrassBlockProps[3];
	        	}

	        	else if (SampledNormals == NORMAL_BOTTOM)
	        	{
	        		data.x = u_GrassBlockProps[7];
	        		data.y = u_GrassBlockProps[8];
	        		data.z = u_GrassBlockProps[9];
	        	}
	        }



            if (u_POM && data.r != u_GrassblockAlbedoID)
            {
                vec2 InitialUV = UV;

                if (SampledNormals == vec3(-1.0f, 0.0f, 0.0f)) {  UV.x = 1.0f - UV.x; UV.y = 1.0f - UV.y; }
                if (SampledNormals == vec3(1.0f, 0.0f, 0.0f)) {  UV.y = 1.0f - UV.y; }
                if (SampledNormals == vec3(0.0f, -1.0f, 0.0f)) {  UV.y = 1.0f - UV.y; }
                //if (SampledNormals == vec3(0.0f, 0.0f, 1.0f)) {  flip = true;}
                //if (SampledNormals == vec3(0.0f, 0.0f, -1.0f)) {  flip = true; }
                
                vec3 TangentViewPosition = tbn * u_ViewerPosition;
                vec3 TangentFragPosition = tbn * WorldPosition.xyz; 
                vec3 TangentViewDirection = normalize(TangentFragPosition - TangentViewPosition);
                
                UV = ParallaxOcclusionMapping(UV, TangentViewDirection, data.z);
                //UV.y = flip && data.r == u_GrassblockAlbedoID ? 1.0f - UV.y : UV.y;
            
            } else { UV.y = 1.0f - UV.y; }

            vec4 PBRMap = texture(u_BlockPBRTextures, vec3(UV, data.z)).rgba;
            vec3 AlbedoColor = texture(u_BlockAlbedoTextures, vec3(UV.xy, data.x)).rgb;
            vec3 NormalMapped = tbn * (texture(u_BlockNormalTextures, vec3(UV, data.y)).rgb * 2.0f - 1.0f);
            vec3 NonAmplifiedNormal = NormalMapped;

            if (u_AmplifyNormalMap) {
                NormalMapped.x *= 1.64f;
                NormalMapped.z *= 1.85f;
                NormalMapped += 1e-4f;
                NormalMapped = normalize(NormalMapped);
            }
            
            float Emissivity = data.w > -0.5f ? textureLod(u_BlockEmissiveTextures, vec3(UV, data.w), 0.0f).r : 0.0f;

            //vec4 Diffuse = BilateralUpsample(u_DiffuseTexture, v_TexCoords, SampledNormals.xyz, WorldPosition.z);
            //vec4 Diffuse = PositionOnlyBilateralUpsample(u_DiffuseTexture, v_TexCoords, WorldPosition.xyz);
            //vec3 Diffuse = DepthOnlyBilateralUpsample(u_DiffuseTexture, v_TexCoords, WorldPosition.z).xyz;
            //vec4 SampledIndirectDiffuse = BilateralUpsample2(u_DiffuseTexture, v_TexCoords, WorldPosition.xyz, SampledNormals.xyz).xyzw;
            
            vec4 SHy;
            vec2 ShCoCg;
            SpatialUpscaleIndirectDiffuse(SampledNormals.xyz, WorldPosition.w, SHy, ShCoCg);

            vec3 IndirectN = NormalMapped.xyz;
            vec3 SampledIndirectDiffuse = SHToIrridiance(SHy, ShCoCg, IndirectN);

            // VXAO from indirect diffuse trace : 
            // 

            bool do_vxao = u_DoVXAO && u_SVGFEnabled;
            if (do_vxao && !u_RTAO && (distance(WorldPosition.xyz, u_ViewerPosition) < 70)) // -> Causes artifacts if the AO is applied too far away
            {
                float bias = 0.00125f;
                if (v_TexCoords.x > bias && v_TexCoords.x < 1.0f - bias &&
                    v_TexCoords.y > bias && v_TexCoords.y < 1.0f - bias)
                {
                    float ao = texture(u_VXAO, v_TexCoords).x;
                    SampledIndirectDiffuse *= vec3(clamp(pow(ao, 1.85f), 0.0f, 1.0f));
                }
            }


            vec3 LightAmbience = (vec3(120.0f, 172.0f, 255.0f) / 255.0f) * 1.01f;
            vec3 Ambient = (AlbedoColor * LightAmbience) * 0.09f;
            float SampledAO = pow(PBRMap.w, 1.25f); // Ambient occlusion map

            // Interpolate and find sun colors : 
            float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
            float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 1.5f), 0.0f, 1.0f);
            vec3 SunColor = mix(SUN_COLOR, DUSK_COLOR * 0.5f, DuskVisibility);

            //vec3 SunColor = SUN_COLOR;
            vec3 SunDirectLighting = CalculateDirectionalLight(WorldPosition.xyz, normalize(u_SunDirection), SunColor, SunColor, AlbedoColor, NormalMapped, PBRMap.xyz, RayTracedShadow);
            vec3 MoonDirectLighting = CalculateDirectionalLight(WorldPosition.xyz, normalize(u_MoonDirection), NIGHT_COLOR, NIGHT_COLOR, AlbedoColor, NormalMapped, PBRMap.xyz, RayTracedShadow);
            vec3 DirectLighting = mix(SunDirectLighting, MoonDirectLighting, SunVisibility * vec3(1.0f));
            
            DirectLighting = (float(!(Emissivity > 0.5f)) * DirectLighting);
            float Roughness = PBRMap.r;
            vec3 DiffuseIndirect = (SampledIndirectDiffuse.xyz * AlbedoColor);
            //vec3 DiffuseIndirect = (SampledIndirectDiffuse.xyz);

            // Compute Specular indirect from spherical harmonic 

            vec4 SpecularSH = texture(u_ReflectionSHData, v_TexCoords);
            vec2 SpecularCoCg = texture(u_ReflectionCoCgData, v_TexCoords).rg;
            vec3 SpecularIndirect = vec3(0.0f);
            
            if (InBiasedSS) {
                SpecularIndirect += SHToIrridiance(SpecularSH, SpecularCoCg, NonAmplifiedNormal); 
            }

            // Dirty hack to make the normals a bit more visible because the reflection map is so low quality 
            // That it hurts my soul
            // (this just slightly darkens areas where the normal map effect is high)
            float NDotMap = pow(abs(dot(SampledNormals.xyz, NormalMapped.xyz)), 200.0f);
            float NDotMapM = clamp(exp(NDotMap) * 0.95f, 0.2f, 1.0f);
            SpecularIndirect *= NDotMapM; 
           

            vec3 Lo = normalize(u_ViewerPosition - WorldPosition.xyz); // Outgoing direction 
            vec3 F0 = mix(vec3(0.04), AlbedoColor, PBRMap.g); // Fresnel at 0 degrees.

            vec3 SpecularFactor = fresnelroughness(Lo, NormalMapped.xyz, vec3(F0), Roughness); 
            
            // Final combine : 
            // Use fresnel to get the amount of diffuse and reflected light 
            DiffuseIndirect *= clamp(SampledAO, 0.2f, 1.01f);

            bool dfg = u_UseDFG;

            if (!dfg) {
                o_Color = (DirectLighting + ((1.0f - SpecularFactor) * DiffuseIndirect) + 
                      (SpecularFactor * SpecularIndirect));
            }

            else {
                float NDotV = clamp(dot(NonAmplifiedNormal, Lo), 0.0f, 1.0f);
                vec3 Tint = mix(vec3(1.0f), AlbedoColor, float(PBRMap.g > 0.0125f));
                vec3 DFG = DFGPolynomialApproximate(F0, Roughness, NDotV); // dfg
                vec3 IndirectSpecularBRDF = DFG * clamp(Tint, 0.0f, 1.0f);

                o_Color = DirectLighting + DiffuseIndirect;
                o_Color += SpecularIndirect * IndirectSpecularBRDF;
            }


             ///if (PBRMap.x > 0.897511) 
            ///{
            ///    o_Color = vec3(1.0f, 0.0f, 0.0f);
            ///}
            ///
            ///else 
            ///{
            ///    o_Color = vec3(0.0f, 1.0f, 0.0f);
            ///}
            
            // Output utility : 
            o_Normal = vec3(NormalMapped.x, NormalMapped.y, NormalMapped.z);
            o_PBR.xyz = PBRMap.xyz;
            o_PBR.w = Emissivity;


            // Fix bloom light leak around the edges : 
            const bool BloomLightLeakFix = true;
            if (BloomLightLeakFix) {
                float lbiasx = 0.032501f;
                float lbiasy = 0.032501f;
                o_PBR.w *= float(UV.x > lbiasx && UV.x < 1.0f - lbiasx &&
                                 UV.y > lbiasy && UV.y < 1.0f - lbiasy);
            }




           

            return;
        }

        else 
        {   
            vec3 CloudAndSky = GetAtmosphereAndClouds();
            o_Color = (CloudAndSky);
            o_Normal = vec3(-1.0f);
            o_PBR.xyz = vec3(-1.0f);
        }
    }

    else 
    {   
        vec3 CloudAndSky = GetAtmosphereAndClouds();
        o_Color = (CloudAndSky);
        o_Normal = vec3(-1.0f);
        o_PBR.xyz = vec3(-1.0f);
    }
}

// cook torrance brdf : 

float ndfGGX(float cosLh, float roughness)
{
	float alpha = roughness * roughness;
	float alphaSq = alpha * alpha;
	float denom = (cosLh * cosLh) * (alphaSq - 1.0) + 1.0;
	return alphaSq / (PI * denom * denom);
}

float gaSchlickG1(float cosTheta, float k)
{
	return cosTheta / (cosTheta * (1.0 - k) + k);
}

float gaSchlickGGX(float cosLi, float cosLo, float roughness)
{
	float r = roughness + 1.0;
	float k = (r * r) / 8.0; // Epic suggests using this roughness remapping for analytic lights.
	return gaSchlickG1(cosLi, k) * gaSchlickG1(cosLo, k);
}

vec3 fresnelSchlick(vec3 F0, float cosTheta)
{
	return F0 + (vec3(1.0) - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 radiance_s, vec3 albedo, vec3 normal, vec3 pbr, float shadow)
{
    const float DIELECTRIC = 0.04f;
    const float Epsilon = 0.00001;
    float Shadow = min(shadow, 1.0f);

    vec3 Lo = normalize(u_ViewerPosition - world_pos);

	vec3 N = normal;
	float cosLo = max(0.0, dot(N, Lo));
	vec3 Lr = 2.0 * cosLo * N - Lo;
	vec3 F0 = mix(vec3(DIELECTRIC), albedo, pbr.g);

    vec3 Li = light_dir;
	vec3 Lradiance = radiance;

	vec3 Lh = normalize(Li + Lo);

	float cosLi = max(0.0, dot(N, Li));
	float cosLh = max(0.0, dot(N, Lh));

	//vec3 F  = fresnelSchlick(F0, max(0.0, dot(Lh, Lo)));
	vec3 F = fresnelroughness(Lo, normal.xyz, vec3(F0), pbr.r); // use fresnel roughness, gives a slightly better approximation
	float D = ndfGGX(cosLh, pbr.r);
	float G = gaSchlickGGX(cosLi, cosLo, pbr.r);

	vec3 kd = mix(vec3(1.0) - F, vec3(0.0), pbr.g);
	vec3 diffuseBRDF = kd * albedo;

	vec3 specularBRDF = (F * D * G) / max(Epsilon, 4.0 * cosLi * cosLo); specularBRDF = clamp(specularBRDF, 0.0f, 2.0f);
	vec3 Result = (diffuseBRDF * Lradiance * cosLi) + (specularBRDF * radiance_s * cosLi);
    return clamp(Result, 0.0f, 2.5) * clamp((1.0f - Shadow), 0.0f, 1.0f);
}

vec4 cubic(float v){
    vec4 n = vec4(1.0, 2.0, 3.0, 4.0) - v;
    vec4 s = n * n * n;
    float x = s.x;
    float y = s.y - 4.0 * s.x;
    float z = s.z - 4.0 * s.y + 6.0 * s.x;
    float w = 6.0 - x - y - z;
    return vec4(x, y, z, w) * (1.0/6.0);
}

vec4 textureBicubic(sampler2D sampler, vec2 texCoords)
{

   vec2 texSize = textureSize(sampler, 0);
   vec2 invTexSize = 1.0 / texSize;

   texCoords = texCoords * texSize - 0.5;


    vec2 fxy = fract(texCoords);
    texCoords -= fxy;

    vec4 xcubic = cubic(fxy.x);
    vec4 ycubic = cubic(fxy.y);

    vec4 c = texCoords.xxyy + vec2 (-0.5, +1.5).xyxy;

    vec4 s = vec4(xcubic.xz + xcubic.yw, ycubic.xz + ycubic.yw);
    vec4 offset = c + vec4 (xcubic.yw, ycubic.yw) / s;

    offset *= invTexSize.xxyy;

    vec4 sample0 = texture(sampler, offset.xz);
    vec4 sample1 = texture(sampler, offset.yz);
    vec4 sample2 = texture(sampler, offset.xw);
    vec4 sample3 = texture(sampler, offset.yw);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix(
       mix(sample3, sample2, sx), mix(sample1, sample0, sx)
    , sy);
}


vec4 TilelessTexture(sampler2DArray samp, vec3 uvt, float v)
{
    vec2 uv = uvt.xy;
    vec2 iuv = floor( uv );
    vec2 fuv = fract( uv );

    // generate per-tile transform
    vec4 ofa = hash4( iuv + vec2(0.0,0.0) );
    vec4 ofb = hash4( iuv + vec2(1.0,0.0) );
    vec4 ofc = hash4( iuv + vec2(0.0,1.0) );
    vec4 ofd = hash4( iuv + vec2(1.0,1.0) );
    
    vec2 ddx = dFdx( uv );
    vec2 ddy = dFdy( uv );

    // transform per-tile uvs
    ofa.zw = sign(ofa.zw-0.5);
    ofb.zw = sign(ofb.zw-0.5);
    ofc.zw = sign(ofc.zw-0.5);
    ofd.zw = sign(ofd.zw-0.5);
    
    // uv's, and derivarives (for correct mipmapping)
    vec2 uva = uv*ofa.zw + ofa.xy; vec2 ddxa = ddx*ofa.zw; vec2 ddya = ddy*ofa.zw;
    vec2 uvb = uv*ofb.zw + ofb.xy; vec2 ddxb = ddx*ofb.zw; vec2 ddyb = ddy*ofb.zw;
    vec2 uvc = uv*ofc.zw + ofc.xy; vec2 ddxc = ddx*ofc.zw; vec2 ddyc = ddy*ofc.zw;
    vec2 uvd = uv*ofd.zw + ofd.xy; vec2 ddxd = ddx*ofd.zw; vec2 ddyd = ddy*ofd.zw;
        
    // fetch and blend
    vec2 b = smoothstep(0.25,0.75,fuv);

    return mix( mix( textureGrad( samp, vec3(uva,uvt.z), ddxa, ddya ), 
                     textureGrad( samp, vec3(uvb,uvt.z), ddxb, ddyb ), b.x ), 
                mix( textureGrad( samp, vec3(uvc,uvt.z), ddxc, ddyc ),
                     textureGrad( samp, vec3(uvd,uvt.z), ddxd, ddyd ), b.x), b.y );
}