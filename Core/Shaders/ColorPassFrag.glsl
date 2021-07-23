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

uniform sampler2D u_DiffuseSHData1;
uniform sampler2D u_DiffuseSHData2;

uniform sampler2D u_ReflectionSHData;
uniform sampler2D u_ReflectionCoCgData;

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

const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

float Noise2d( in vec2 x )
{
    float xhash = cos( x.x * 37.0 );
    float yhash = cos( x.y * 57.0 );
    return fract( 415.92653 * ( xhash + yhash ) );
}

float NoisyStarField( in vec2 vSamplePos, float fThreshhold )
{
    float StarVal = Noise2d( vSamplePos );
    if ( StarVal >= fThreshhold )
        StarVal = pow( (StarVal - fThreshhold)/(1.0 - fThreshhold), 6.0 );
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

// Original star shader by : https://www.shadertoy.com/view/Md2SR3
float StableStarField( in vec2 vSamplePos, float fThreshhold )
{
    float fractX = fract( vSamplePos.x );
    float fractY = fract( vSamplePos.y );
    vec2 floorSample = floor( vSamplePos );
    float v1 = NoisyStarField( floorSample, fThreshhold );
    float v2 = NoisyStarField( floorSample + vec2( 0.0, 1.0 ), fThreshhold );
    float v3 = NoisyStarField( floorSample + vec2( 1.0, 0.0 ), fThreshhold );
    float v4 = NoisyStarField( floorSample + vec2( 1.0, 1.0 ), fThreshhold );

    float StarVal =   v1 * ( 1.0 - fractX ) * ( 1.0 - fractY )
        			+ v2 * ( 1.0 - fractX ) * fractY
        			+ v3 * fractX * ( 1.0 - fractY )
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

    float star = StableStarField(uv * 700.0f, 0.999);
    star *= 2.0f;
    float rand_val = rand(fragpos.xy);
	return clamp(star, 0.0f, 100000.0f) * 30.0f;
}

vec2 RayBoxIntersect(vec3 boundsMin, vec3 boundsMax, vec3 rayOrigin, vec3 invRaydir)
{
	vec3 t0 = (boundsMin - rayOrigin) * invRaydir;
	vec3 t1 = (boundsMax - rayOrigin) * invRaydir;
	vec3 tmin = min(t0, t1);
	vec3 tmax = max(t0, t1);
	
	float dstA = max(max(tmin.x, tmin.y), tmin.z);
	float dstB = min(tmax.x, min(tmax.y, tmax.z));
	
	// CASE 1: ray intersects box from outside (0 <= dstA <= dstB)
	// dstA is dst to nearest intersection, dstB dst to far intersection
	
	// CASE 2: ray intersects box from inside (dstA < 0 < dstB) 
	// dstA is the dst to intersection behind the ray, dstB is dst to forward intersection
	
	// CASE 3: ray misses box (dstA > dstB)
	
	float dstToBox = max(0, dstA);
	float dstInsideBox = max(0, dstB - dstToBox);
	return vec2(dstToBox, dstInsideBox);
}

bool GetAtmosphere(inout vec3 atmosphere_color, in vec3 in_ray_dir)
{
    vec3 sun_dir = normalize(u_SunDirection); 
    vec3 moon_dir = vec3(-sun_dir.x, -sun_dir.y, sun_dir.z); 

    vec3 ray_dir = normalize(in_ray_dir);
    
    if(dot(ray_dir, sun_dir) > 0.999825f)
    {
        atmosphere_color = ATMOSPHERE_SUN_COLOR; return true;
    }

    if(dot(ray_dir, moon_dir) > 0.99986f)
    {
        atmosphere_color = ATMOSPHERE_MOON_COLOR; return true;
    }

    vec3 atmosphere = texture(u_Skybox, ray_dir).rgb;

    float star_visibility;
    star_visibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; 
    star_visibility = 1.0f - star_visibility;
    vec3 stars = vec3(stars(vec3(in_ray_dir)) * star_visibility);
    stars = clamp(stars, 0.0f, 1.3f);

    atmosphere += stars;

    atmosphere_color = atmosphere;

    return false;
}

vec3 GetAtmosphereAndClouds(vec3 Sky)
{
    if (!u_CloudsEnabled)
    {
        return Sky;
    }

    float BoxSize = (u_CloudBoxSize - 8.0f); // -8 to reduce artifacts.
    vec3 origin = vec3(v_RayOrigin.x, 0.0f, v_RayOrigin.z);
    vec2 Dist = RayBoxIntersect(origin + vec3(-BoxSize, CLOUD_HEIGHT, -BoxSize), origin + vec3(BoxSize, CLOUD_HEIGHT - 12, BoxSize), v_RayOrigin, 1.0f / (v_RayDirection));
    bool Intersect = !(Dist.y == 0.0f);

    if (!Intersect)
    {
        return Sky;
    }

    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
    const vec3 D = (vec3(355.0f, 10.0f, 0.0f) / 255.0f) * 0.4f;
    vec3 S = vec3(1.45f);
    float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 1.8f), 0.0f, 1.0f);
    S = mix(S, D, DuskVisibility);
    vec3 M = mix(S + 0.001f, (vec3(46.0f, 142.0f, 255.0f) / 255.0f) * 0.1f, SunVisibility); 
	vec3 IntersectionPosition = v_RayOrigin + (normalize(v_RayDirection) * Dist.x);
	vec4 SampledCloudData = texture(u_CloudData, v_TexCoords).rgba;
    vec3 Scatter = SampledCloudData.xyz;
    float Transmittance = SampledCloudData.w;
    return (Sky * 1.0f) * clamp(Transmittance, 0.95f, 1.0f) + (Scatter * 1.0f * M); // see ya pbr
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

    vec3 BiasedWorldPos = world_pos - (flat_normal * 0.125f);
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

	if (normal == Normals[0])
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[0];
		bitangent = BiTangents[0];
    }

    else if (normal == Normals[1])
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[1];
		bitangent = BiTangents[1];
    }

    else if (normal == Normals[2])
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[2];
		bitangent = BiTangents[2];
    }

    else if (normal == Normals[3])
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[3];
		bitangent = BiTangents[3];
    }
	
    else if (normal == Normals[4])
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[4];
		bitangent = BiTangents[4];
    }
    

    else if (normal == Normals[5])
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[5];
		bitangent = BiTangents[5];
    }
}

vec3 fresnelroughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
{
	return F0 + (max(vec3(pow(1.0f - roughness, 3.0f)) - F0, vec3(0.0f))) * pow(max(1.0 - clamp(dot(Eye, norm), 0.0f, 1.0f), 0.0f), 5.0f);
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

// COLORS //
const vec3 SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * 4.0f;
const vec3 NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.25f; 
const vec3 DUSK_COLOR = (vec3(255.0f, 204.0f, 144.0f) / 255.0f) * 0.064f; 

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

    vec3 SampledNormals = texture(u_NormalTexture, v_TexCoords).rgb;
    vec3 AtmosphereAt = vec3(0.0f);

    o_Color = vec3(1.0f);
    bool BodyIntersect = GetAtmosphere(AtmosphereAt, v_RayDirection);
    o_PBR.w = float(BodyIntersect);

    if (WorldPosition.w > 0.0f)
    {
        float RayTracedShadow = 0.0f;
        
        RayTracedShadow = ComputeShadow(WorldPosition.xyz, SampledNormals.xyz);

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
            
            if (u_AmplifyNormalMap) {
                NormalMapped.x *= 1.64f;
                NormalMapped.z *= 1.85f;
                NormalMapped += 1e-4f;
                NormalMapped = normalize(NormalMapped);
            }
            
            float Emissivity = data.w > -0.5f ? texture(u_BlockEmissiveTextures, vec3(UV, data.w)).r : 0.0f;

            //vec4 Diffuse = BilateralUpsample(u_DiffuseTexture, v_TexCoords, SampledNormals.xyz, WorldPosition.z);
            //vec4 Diffuse = PositionOnlyBilateralUpsample(u_DiffuseTexture, v_TexCoords, WorldPosition.xyz);
            //vec3 Diffuse = DepthOnlyBilateralUpsample(u_DiffuseTexture, v_TexCoords, WorldPosition.z).xyz;
            //vec4 SampledIndirectDiffuse = BilateralUpsample2(u_DiffuseTexture, v_TexCoords, WorldPosition.xyz, SampledNormals.xyz).xyzw;
            
            vec4 SHy = texture(u_DiffuseSHData1, v_TexCoords);
            vec2 ShCoCg = texture(u_DiffuseSHData2, v_TexCoords).xy;

            vec3 IndirectN = NormalMapped.xyz;
            vec3 SampledIndirectDiffuse = vec3(0.0f);
            int SampleCount = 4;

            for (int i = 0 ; i < SampleCount ; i++) {
                vec3 SampleDirection = IndirectN;
                vec2 J = vec2(nextFloat(RNG_SEED, 1.0f, 1.2f), nextFloat(RNG_SEED, 1.0f, 1.165f));
                SampleDirection.x *= J.x;
                SampleDirection.z *= J.y;
                SampleDirection = normalize(SampleDirection);
                SampledIndirectDiffuse += SHToIrridiance(SHy, ShCoCg, SampleDirection);
            }
            
            SampledIndirectDiffuse /= float(SampleCount);
            //vec3 SampledIndirectDiffuse = SHToIrridiance(SHy, ShCoCg, IndirectN);


            //float AO = texture(u_DiffuseTexture, v_TexCoords).w;
            float AO = 1.0f;

            vec3 LightAmbience = (vec3(120.0f, 172.0f, 255.0f) / 255.0f) * 1.01f;
            vec3 Ambient = (AlbedoColor * LightAmbience) * 0.09f;
            float SampledAO = pow(PBRMap.w, 1.25f); // Ambient occlusion map

            float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
            float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 1.5f), 0.0f, 1.0f);
            vec3 SunColor = mix(SUN_COLOR, DUSK_COLOR * 0.5f, DuskVisibility);
            //vec3 SunColor = SUN_COLOR;
            vec3 SunDirectLighting = CalculateDirectionalLight(WorldPosition.xyz, normalize(u_SunDirection), SunColor, SunColor * 0.4f, AlbedoColor, NormalMapped, PBRMap.xyz, RayTracedShadow);
            vec3 MoonDirectLighting = CalculateDirectionalLight(WorldPosition.xyz, normalize(u_MoonDirection), NIGHT_COLOR, NIGHT_COLOR, AlbedoColor, NormalMapped, PBRMap.xyz, RayTracedShadow);
            vec3 DirectLighting = mix(SunDirectLighting, MoonDirectLighting, SunVisibility * vec3(1.0f));
            
            DirectLighting = (float(!(Emissivity > 0.5f)) * DirectLighting);
            float Roughness = PBRMap.r;
            vec3 DiffuseIndirect = (SampledIndirectDiffuse.xyz * AlbedoColor);

            // Compute Specular indirect from spherical harmonic 

            vec4 SpecularSH = texture(u_ReflectionSHData, v_TexCoords);
            vec2 SpecularCoCg = texture(u_ReflectionCoCgData, v_TexCoords).rg;
            vec3 SpecularIndirect = SHToIrridiance(SpecularSH, SpecularCoCg, IndirectN).rgb;
            
            // Dirty hack to make the normals a bit more visible because the reflection map is so low quality 
            // That it hurts my soul
            float NDotMap = pow(abs(dot(SampledNormals.xyz, NormalMapped.xyz)), 200.0f);
            float NDotMapM = clamp(exp(NDotMap) * 0.95f, 0.2f, 1.0f);
            SpecularIndirect *= NDotMapM; 
           
            vec3 Lo = normalize(u_ViewerPosition - WorldPosition.xyz); // Outgoing direction 
            vec3 F0 = mix(vec3(0.04), AlbedoColor, PBRMap.g); // Fresnel at 0 degrees.

            vec3 SpecularFactor = fresnelroughness(Lo, NormalMapped.xyz, vec3(F0), Roughness); 
            
            o_Color = (DirectLighting + ((1.0f - SpecularFactor) * DiffuseIndirect) + 
                      (SpecularFactor * SpecularIndirect * min((PBRMap.g + 1.0f), 1.3f))) 
                      * clamp(SampledAO, 0.2f, 1.01f);

            o_Normal = vec3(NormalMapped.x, NormalMapped.y, NormalMapped.z);
            o_PBR.xyz = PBRMap.xyz;
            o_PBR.w = Emissivity;

            const bool BloomLightLeakFix = true;
            if (BloomLightLeakFix) {
                float lbiasx = 0.02501f;
                float lbiasy = 0.03001f;
                o_PBR.w *= float(UV.x > lbiasx && UV.x < 1.0f - lbiasx &&
                                 UV.y > lbiasy && UV.y < 1.0f - lbiasy);
            }

            if (!u_RTAO && (distance(WorldPosition.xyz, u_ViewerPosition) < 80)) // -> Causes artifacts if the AO is applied too far away
            {
                float bias = 0.00125f;
                if (v_TexCoords.x > bias && v_TexCoords.x < 1.0f - bias &&
                    v_TexCoords.y > bias && v_TexCoords.y < 1.0f - bias)
                {
                    o_Color *= vec3(clamp(pow(AO, 1.0f), 0.0f, 1.0f));
                }
            }

            return;
        }

        else 
        {   
            vec3 CloudAndSky = GetAtmosphereAndClouds(AtmosphereAt);
            o_Color = (CloudAndSky);
            o_Normal = vec3(-1.0f);
            o_PBR.xyz = vec3(-1.0f);
            o_PBR.w = float(BodyIntersect) * 0.35f;
        }
    }

    else 
    {   
        vec3 CloudAndSky = GetAtmosphereAndClouds(AtmosphereAt);
        o_Color = (CloudAndSky);
        o_Normal = vec3(-1.0f);
        o_PBR.xyz = vec3(-1.0f);
        o_PBR.w = float(BodyIntersect) * 0.35f;
    }
}

float ndfGGX(float cosLh, float roughness)
{
	float alpha   = roughness * roughness;
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
    const float Epsilon = 0.00001;
    float Shadow = min(shadow, 1.0f);

    vec3 Lo = normalize(u_ViewerPosition - world_pos);

	vec3 N = normal;
	float cosLo = max(0.0, dot(N, Lo));
	vec3 Lr = 2.0 * cosLo * N - Lo;
	vec3 F0 = mix(vec3(0.04), albedo, pbr.g);

    vec3 Li = light_dir;
	vec3 Lradiance = radiance;

	vec3 Lh = normalize(Li + Lo);

	float cosLi = max(0.0, dot(N, Li));
	float cosLh = max(0.0, dot(N, Lh));

	//vec3 F  = fresnelSchlick(F0, max(0.0, dot(Lh, Lo)));
	vec3 F = fresnelroughness(Lo, normal.xyz, vec3(F0), pbr.r); 
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