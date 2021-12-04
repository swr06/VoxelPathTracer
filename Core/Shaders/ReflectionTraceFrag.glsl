#version 430 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

//#define DERIVE_FROM_DIFFUSE_SH

#define REPROJECT_TO_SCREEN_SPACE
#define TRACE_LENGTH 64
#define clamp01(x) (clamp(x,0.0f,1.0F))

//#define ALBEDO_TEX_LOD 3 // 512, 256, 128
//#define JITTER_BASED_ON_ROUGHNESS

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


layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out float o_HitDistance;
layout (location = 3) out float o_EmissivityHitMask; // -> Used as input to the firefly rejection filter 

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_InitialTraceNormalTexture;
uniform sampler2D u_BlockIDTex;
uniform sampler2D u_DataTexture;

uniform float u_SunStrengthModifier;
uniform float u_MoonStrengthModifier;

uniform sampler2D u_ProjectedClouds;

//uniform sampler2D u_BlueNoiseTexture;

uniform sampler3D u_VoxelData;
uniform sampler3D u_DistanceFieldTexture;

uniform sampler2D u_PlayerSprite;

uniform sampler2D u_DiffuseSH;
uniform sampler2D u_DiffuseCoCg;
uniform sampler2D u_ShadowTrace;

uniform sampler2D u_IndirectAO;

uniform sampler3D u_LPV;
uniform usampler3D u_LPVBlocks;

uniform bool u_CloudReflections;


uniform bool u_TemporalFilterReflections;


uniform bool u_RoughReflections;
uniform bool CHECKERBOARD_SPEC_SPP;

uniform samplerCube u_Skymap;

uniform float u_ReflectionTraceRes;

uniform vec2 u_Dimensions;
uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;
uniform vec3 u_StrongerLightDirection;
uniform float u_Time;

uniform int u_GrassBlockProps[10];

uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlockEmissiveTextures;

uniform vec3 u_ViewerPosition;
uniform int u_SPP;

uniform int u_CurrentFrame;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform bool u_ReprojectToScreenSpace;
uniform bool u_UseBlueNoise;

uniform int u_CurrentFrameMod128;


uniform bool TEMPORAL_SPEC=false;
uniform bool u_ReflectPlayer;

// LPVGI
uniform bool u_LPVGI;
uniform bool u_QualityLPVGI;

uniform vec2 u_Halton;
uniform bool u_RoughnessBias;


layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

layout (std430, binding = 2) buffer BlueNoise_Data
{
	int sobol_256spp_256d[256*256];
	int scramblingTile[128*128*8];
	int rankingTile[128*128*8];
};


layout (std430, binding = 4) buffer SSBO_BlockAverageData
{
    vec4 BlockAverageColorData[128]; // Returns the average color per block type 
};




const vec3 ATMOSPHERE_SUN_COLOR = vec3(1.0f * 6.25f, 1.0f * 6.25f, 0.8f * 4.0f);
const vec3 ATMOSPHERE_MOON_COLOR =  vec3(0.7f, 0.7f, 1.25f);
const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);
vec3 SkyAmbientG = vec3(0.0f);
vec3 SAMPLED_SUN_COLOR = vec3(0.0f);
vec3 SAMPLED_MOON_COLOR = vec3(0.0f);
vec3 SAMPLED_COLOR_MIXED = vec3(0.0f);
vec3 C_SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * 6.0f * u_SunStrengthModifier;
vec3 C_NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * vec3(0.9,0.9,1.0f) * 0.225f * u_MoonStrengthModifier; 
vec3 C_DUSK_COLOR = (vec3(255.0f, 204.0f, 144.0f) / 255.0f) * 0.1f; 

		
// Function prototypes
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType, bool shadow);
float GetVoxel(ivec3 loc);
float GetShadowAt(in vec3 pos, in vec3 ldir);
void ComputePlayerReflection(in vec3 ro, in vec3 rd, inout vec3 col, float block_t);


float Luma(vec3 x) { return dot(x, vec3(0.2125, 0.7154, 0.0721)); }



float samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp(ivec2 px, int sampleIndex, int sampleDimension)
{
	int pixel_i = px.x;
	int pixel_j = px.y;
	
	// wrap arguments
	pixel_i = pixel_i & 127;
	pixel_j = pixel_j & 127;
	sampleIndex = sampleIndex & 255;
	sampleDimension = sampleDimension & 255;

	// xor index based on optimized ranking
	int rankedSampleIndex = sampleIndex ^ rankingTile[sampleDimension + (pixel_i + pixel_j*128)*8];

	// fetch value in sequence
	int value = sobol_256spp_256d[sampleDimension + rankedSampleIndex*256];

	// If the dimension is optimized, xor sequence value based on optimized scrambling
	value = value ^ scramblingTile[(sampleDimension%8) + (pixel_i + pixel_j*128)*8];

	// convert to float and return
	float v = (0.5f+value)/256.0f;
	return v;
}

vec3 BasicSaturation(vec3 Color, float Adjustment)
{
    const vec3 LuminosityCoefficients = vec3(0.2125f, 0.7154f, 0.0721f);
    vec3 Luminosity = vec3(dot(Color, LuminosityCoefficients));
    return mix(Luminosity, Color, Adjustment);
}

// GGX ->

float ndfGGX(float cosLh, float roughness)
{
	float alpha   = roughness * roughness;
	float alphaSq = alpha * alpha;

	float denom = (cosLh * cosLh) * (alphaSq - 1.0) + 1.0;
	return alphaSq / (PI * denom * denom);
}

// Fresnel ->

float gaSchlickG1(float cosTheta, float k)
{
	return cosTheta / (cosTheta * (1.0 - k) + k);
}

// Geometry ->

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

// Cook torrance brdf ->
vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 albedo, vec3 normal, vec3 pbr, float shadow)
{
    const float Epsilon = 0.00001;
    float Shadow = min(shadow, 1.0f);
    vec3 Lo = normalize(u_ViewerPosition - world_pos);
	vec3 N = normal;
	float cosLo = max(0.0, dot(N, Lo));
	vec3 Lr = 2.0 * cosLo * N - Lo;
	vec3 F0 = mix(vec3(0.04), albedo, pbr.g);
    vec3 Li = light_dir; // no need to normalize
	vec3 Lradiance = radiance;
	vec3 Lh = normalize(Li + Lo);
	float cosLi = max(0.0, dot(N, Li));
	float cosLh = max(0.0, dot(N, Lh));
	vec3 F  = fresnelSchlick(F0, max(0.0, dot(Lh, Lo)));
	float D = ndfGGX(cosLh, pbr.r);
	float G = gaSchlickGGX(cosLi, cosLo, pbr.r);
	vec3 kd = mix(vec3(1.0) - F, vec3(0.0), pbr.g);
	vec3 diffuseBRDF = kd * albedo;
	vec3 specularBRDF = (F * D * G) / max(Epsilon, 4.0 * cosLi * cosLo);
	vec3 radiance_s = radiance * 0.05f * 0;
	vec3 Result = (diffuseBRDF * Lradiance * cosLi) + (specularBRDF * radiance_s * cosLi);
    return max(Result, 0.0f) * clamp((1.0f - Shadow), 0.0f, 1.0f);
}

// Samples the ggx vndf and returns a microfacet normal ->
vec3 ImportanceSampleGGX(vec3 N, float roughness, vec2 Xi)
{
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
	
    float phi = 2.0 * PI * Xi.x;
    float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (alpha2 - 1.0) * Xi.y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	
    vec3 H;
    H.x = cos(phi) * sinTheta;
    H.y = sin(phi) * sinTheta;
    H.z = cosTheta;
	
    vec3 up = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(up, N));
    vec3 bitangent = cross(N, tangent);
	
    vec3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
    return normalize(sampleVec);
} 


// By iq
vec4 TextureSmooth(sampler2D samp, vec2 uv) 
{
    vec2 textureResolution = textureSize(samp, 0).xy;
	uv = uv*textureResolution + 0.5f;
	vec2 iuv = floor(uv);
	vec2 fuv = fract(uv);
	uv = iuv + fuv*fuv*(3.0f-2.0f*fuv); 
	uv = (uv - 0.5f) / textureResolution;
	return texture(samp, uv).xyzw;
}

vec2 ProjectDirection(vec3 Direction, vec2 TextureSize);

vec3 RetrieveProjectedClouds(vec3 Sky, vec3 R)
{
	if (!u_CloudReflections) {
		return Sky;
	}

	vec2 TextureSize = vec2(textureSize(u_ProjectedClouds,0).xy);
	vec2 ProjectedDirection = ProjectDirection(R, TextureSize);

	if (ProjectedDirection != clamp(ProjectedDirection, 0.001f, 0.999f)) {
		return Sky;
	}
	
	// Compute some fake colors and look close enough ->
	float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
    vec3 BaseSun = SAMPLED_SUN_COLOR * vec3(1.0f, 0.9f, 0.9f) * 1.0f;
	BaseSun *= 0.6f;
    float DuskVisibility = clamp(pow(abs(u_SunDirection.y - 1.0), 2.0f), 0.0f, 1.0f);
    BaseSun = mix(vec3(1.5f, 1.5f, 1.55f), BaseSun, DuskVisibility);
	BaseSun = clamp01(BaseSun);
    vec3 ScatterColor = mix(BaseSun, BasicSaturation(SAMPLED_MOON_COLOR, 0.5f) * 1.32525f, SunVisibility); 
	ScatterColor = clamp01(ScatterColor);
	ScatterColor *= 0.8f;

	// Fake ambient ->
	vec3 SkyAmbient = pow(max(SkyAmbientG, 0.45f), vec3(1.25f)) * 1.75f;

	// Fetch with smooth filter to reduce bilinear artifacts ->
	vec4 CloudFetch = TextureSmooth(u_ProjectedClouds, ProjectedDirection);

	// Use projected clouds
	vec3 Return = Sky * max(0.4f, CloudFetch.w) + clamp((CloudFetch.xyz * ScatterColor * vec3(0.8f, 0.8f, 1.0f) * SkyAmbient * 1.2f), 0.0f, 1.0f);
	return Return;
}	


void GetAtmosphere(inout vec3 Out, in vec3 V)
{
    vec3 SunDirection = (u_SunDirection); 
    vec3 MoonDirection = vec3(-SunDirection.x, -SunDirection.y, SunDirection.z); 
    vec3 NormalizedV = normalize(V);
    vec3 SampledSky = texture(u_Skymap, NormalizedV).rgb;
	SampledSky = RetrieveProjectedClouds(SampledSky, NormalizedV);
    Out = SampledSky;
}


bool GetPlayerIntersect(in vec3 pos, in vec3 ldir);

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

// Quake2RTX
float[6] IrridianceToSH(vec3 Radiance, vec3 Direction) {
	
	float Co = Radiance.x - Radiance.z; 
	float T = Radiance.z + Co * 0.5f; 
	float Cg = Radiance.y - T;
	float Y  = max(T + Cg * 0.5f, 0.0);
	float L00  = 0.282095f;
    float L1_1 = 0.488603f * Direction.y;
    float L10  = 0.488603f * Direction.z;
    float L11  = 0.488603f * Direction.x;
	float ReturnValue[6];
	ReturnValue[0] = max(L11 * Y, -100.0f);
	ReturnValue[1] = max(L1_1 * Y, -100.0f);
	ReturnValue[2] = max(L10 * Y, -100.0f);
	ReturnValue[3] = max(L00 * Y, -100.0f);

	ReturnValue[4] = Co;
	ReturnValue[5] = Cg;
	return ReturnValue;
}

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

vec3 SHToIrradianceA(vec4 shY, vec2 CoCg)
{
	float Y = max(0, 3.544905f * shY.w);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

// Shift x texture coordinate if the current step is a checker step ->
vec2 GetCheckerboardedUV()
{
	vec2 Screenspace = v_TexCoords;
	float TexelSizeX = 1.0f / u_Dimensions.x;
	Screenspace.x += (float(int(gl_FragCoord.x + gl_FragCoord.y) % 2 == int(u_CurrentFrame % 2)) * TexelSizeX);
	return Screenspace;
}

// Fetch correct texture ids 
vec4 GetTextureIDs(int BlockID) 
{
	return vec4(float(BlockAlbedoData[BlockID]),
				float(BlockNormalData[BlockID]),
				float(BlockPBRData[BlockID]),
				float(BlockEmissiveData[BlockID]));
}

// Samples the block id data texture 
int GetBlockID(vec2 txc)
{
	float id = texture(u_BlockIDTex, txc).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

bool InThresholdedScreenSpace(in vec2 v) 
{
	float b = 0.032593f;
	return v.x > b && v.x < 1.0f - b && v.y > b && v.y < 1.0f - b;
}

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

vec3 SampleNormalFromTex(sampler2D samp, vec2 txc) { 
    return GetNormalFromID(texture(samp, txc).x);
}



float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}



// Reprojects hit coordinate to screen space to try and infer indirect lighting and the visibility term for the direct brdf ->
vec2 ReprojectReflectionToScreenSpace(vec3 HitPosition, vec3 HitNormal, out bool Success, out vec4 PositionAt)
{
    vec4 ProjectedPosition = u_Projection * u_View * vec4(HitPosition, 1.0f);
    ProjectedPosition.xyz /= ProjectedPosition.w;
    ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;
    PositionAt = GetPositionAt(u_PositionTexture, ProjectedPosition.xy);
    vec3 NormalAt = SampleNormalFromTex(u_InitialTraceNormalTexture, ProjectedPosition.xy).xyz;
    vec3 PositionDifference = abs(PositionAt.xyz - HitPosition.xyz);
    float Error = dot(PositionDifference, PositionDifference);
    Success = Error < 0.095f && NormalAt == HitNormal && InThresholdedScreenSpace(ProjectedPosition.xy);
    return ProjectedPosition.xy;
}

// rng ->
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

// Samples blue noise data ->
float SampleBlueNoise1D(ivec2 Pixel, int Index, int Dimension) {
	return samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp(Pixel, Index, Dimension);
}

int CurrentBLSample = 0;

vec2 SampleBlueNoise2D(int Index) 
{
	vec2 Noise;
	Noise.x = SampleBlueNoise1D(ivec2(gl_FragCoord.xy), Index, 1+CurrentBLSample);
	Noise.y = SampleBlueNoise1D(ivec2(gl_FragCoord.xy), Index, 2+CurrentBLSample);
	CurrentBLSample += 2;
	CurrentBLSample = CurrentBLSample % 128;
	return Noise;
}


// Gets the most closest reflection direction to help reduce variance 
// *NOT* physically correct as this will add bias which is not accounted for in the pdf 
vec3 GetReflectionDirection(vec3 N, float R) {

	// Gets a reflection vector that is closest to the normal
	R = max(R, 0.05f);
	float NearestDot = -100.0f;
	vec3 BestDirection;

	for (int i = 0 ; i < 4 ; i++) {
		vec2 Xi = u_UseBlueNoise ? SampleBlueNoise2D(TEMPORAL_SPEC?u_CurrentFrameMod128:100) : hash2();
		Xi = Xi * vec2(1.0f, 0.75f);
		vec3 ImportanceSampled = ImportanceSampleGGX(N, R, Xi);
		float d = dot(ImportanceSampled,N);
		if (d > NearestDot) {
			BestDirection = ImportanceSampled;
			NearestDot = d;
		}
	}

	return BestDirection;
}

bool SunStronger;

vec3 TemperatureToRGB(float temperatureInKelvins);

vec3 SampleSunColor()
{
    const vec3 TemperatureModifier = TemperatureToRGB(5778.0f);
    vec3 SunTransmittance = texture(u_Skymap, u_SunDirection.xyz).xyz;
    vec3 SunColor = SunTransmittance;
    SunColor *= TemperatureModifier;
    return SunColor * PI * 2.2f * u_SunStrengthModifier;
}

vec3 SampleMoonColor()
{
    vec3 MoonTransmittance = texture(u_Skymap, u_MoonDirection).xyz;
    vec3 MoonColor = MoonTransmittance;
    MoonColor = MoonColor * PI * u_MoonStrengthModifier;
    float L = GetLuminance(MoonColor);
    MoonColor = BasicSaturation(MoonColor, 1.3f); 
    return MoonColor * 0.42525f * u_MoonStrengthModifier;
}

vec3 SampleLPVData(vec3 UV);


// Samples LPV for GI
vec3 SampleLPV(vec3 P, vec3 B)
{
	vec3 Sky = SkyAmbientG; 
	float L = dot(Sky, vec3(0.2125, 0.7154, 0.0721));

	// Fake it till you make it sky gi

	// Fake desat ->
	Sky = mix(vec3(L), Sky, SunStronger ? 0.3f : 0.6f); 
	vec3 Skylighting = texture(u_IndirectAO, v_TexCoords).y * (SunStronger ? 3.5f : 4.0f) * Sky; // The y component of the ao texture contains the amount of skylight
	Skylighting += bayer16(gl_FragCoord.xy)/128.0f;
	Skylighting = clamp(Skylighting, vec3(0.0f), B + vec3(0.075f));

	vec3 LPV = SampleLPVData(P);
	return Skylighting + LPV;
}

vec3 LPVDither = vec3(0.0f);


void main()
{
	SkyAmbientG = texture(u_Skymap, vec3(0.0f, 1.0f, 0.0f)).xyz;

	LPVDither = vec3(bayer32(gl_FragCoord.xy+vec2(u_CurrentFrame*0.75,u_CurrentFrame*0.5)*float(u_TemporalFilterReflections))); //texture(u_BlueNoiseHighRes, (g_TexCoords * 0.5f * (u_Dimensions / vec2(textureSize(u_BlueNoiseHighRes,0).xy)))).xyz;
    const vec3 VolumeResolution = vec3(384.0f, 128.0f, 384.0f);
    LPVDither /= VolumeResolution;


	SAMPLED_SUN_COLOR = SampleSunColor();
    SAMPLED_MOON_COLOR = SampleMoonColor();
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
    SAMPLED_COLOR_MIXED = mix(SAMPLED_SUN_COLOR, SAMPLED_MOON_COLOR, SunVisibility);


	float TIME = TEMPORAL_SPEC ? u_Time : 1.0f;
	//RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * 800 * int(floor(fract(TIME) * 200));
	//RNG_SEED ^= RNG_SEED << 13;
    //RNG_SEED ^= RNG_SEED >> 17;
    //RNG_SEED ^= RNG_SEED << 5;
	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(TIME) * 100.0f;

	// checker spp
	bool CheckerStep = int(gl_FragCoord.x + gl_FragCoord.y) % 2 == (u_CurrentFrame % 2);
	int SPP = clamp(u_SPP, 1, 16);
	if (CHECKERBOARD_SPEC_SPP) {
		SPP = int(mix(u_SPP, (u_SPP + u_SPP % 2) / 2, float(CheckerStep)));
	}

	SPP = clamp(SPP, 1, 16);

	int total_hits = 0;

	vec4 TotalSH = vec4(0.0f);
	vec2 TotalCoCg = vec2(0.0f);

	// Jitter ->
	vec2 g_TexCoords = v_TexCoords;
	vec2 Jitter = u_Halton;
	Jitter = clamp(Jitter * 1.0f, -2.0f, 2.0f);
	vec2 JitteredUV = g_TexCoords + ((Jitter / u_Dimensions) * float(u_TemporalFilterReflections));
	//float JitteredGBuffer = texelFetch(u_PositionTexture, ivec2(JitteredUV * textureSize(u_PositionTexture, 0)), 0).x;
	if (true) {
		g_TexCoords = JitteredUV;
	}



	// Current gbuffer texel is the sky ->
	vec4 SampledWorldPosition = GetPositionAt(u_PositionTexture, g_TexCoords); // initial intersection point
	
	if (SampledWorldPosition.w < 0.0f)
	{
		o_SH = vec4(0.0f);
		o_CoCg.xy = vec2(0.0f, 0.0f);
		o_EmissivityHitMask = 0.0f;
		o_HitDistance = -1.0f;
		return;
	}


	vec3 InitialTraceNormal = SampleNormalFromTex(u_InitialTraceNormalTexture, g_TexCoords).rgb;
	vec4 data = GetTextureIDs(GetBlockID(g_TexCoords));

	vec2 iUV; 
	vec3 iTan, iBitan;
	CalculateVectors(SampledWorldPosition.xyz, InitialTraceNormal, iTan, iBitan, iUV);
	iUV.y = 1.0f - iUV.y;
	vec4 PBRMap = texture(u_BlockPBRTextures, vec3(iUV, data.z)).rgba; // -> Base pbr data 


	float RoughnessAt = PBRMap.r;

	float MetalnessAt = PBRMap.g;
	vec3 I = normalize(SampledWorldPosition.xyz - u_ViewerPosition); // Incident 
	mat3 tbn = mat3((iTan), (iBitan), (InitialTraceNormal)); 
	
	vec3 NormalMappedInitial = tbn * (texture(u_BlockNormalTextures, vec3(vec2(iUV.x, iUV.y), data.g)).rgb * 2.0f - 1.0f);
	SampledWorldPosition.xyz += InitialTraceNormal.xyz * 0.04f; // Apply bias.
    NormalMappedInitial = normalize(NormalMappedInitial);
	
	float ComputedShadow = 0.0f;
	int ShadowItr = 0;

	vec3 NormalizedStrongerDir = (u_StrongerLightDirection);

	float MaxHitDistance = -1.0f;
	bool Hit = false;
	vec3 ReflectionVector;

	//int MaxSPP = SPP;
	//int MinSPP = clamp(SPP / 2, 2, 32);
	//SPP = 4;

	vec3 refpos = SampledWorldPosition.xyz - (InitialTraceNormal * 0.5f);  // Bias 

	vec4 DiffuseSH = texture(u_DiffuseSH, v_TexCoords).rgba;
	vec2 DiffuseCoCg = texture(u_DiffuseCoCg, v_TexCoords).rg;
	vec3 BaseIndirectDiffuse = SHToIrradianceA(DiffuseSH, DiffuseCoCg);

	float AveragedHitDistance = 0.001f;
	float TotalMeaningfulHits = 0.0f;
	
	SunStronger = u_StrongerLightDirection == u_SunDirection;
	

	float EmissivityMask = 0.0f;

	// Roughness bias ->
	float RoughnessBias = 1.0f;
	RoughnessBias = mix(1.0f, 0.85, float(u_RoughnessBias));


	for (int s = 0 ; s < SPP ; s++)
	{
		//if (MetalnessAt < 0.025f) 
		//{
		//	continue;
		//}

		// importance sample :
		vec3 ReflectionNormal = u_RoughReflections ? GetReflectionDirection(NormalMappedInitial,clamp(RoughnessAt*RoughnessBias,0.01,1.)) : NormalMappedInitial;
		vec3 R = (reflect(I, ReflectionNormal)); ReflectionVector = R;

		#ifdef DERIVE_FROM_DIFFUSE_SH
		if (RoughnessAt >= 0.85f) { 
			vec3 DiffuseSHCompute = SHToIrridiance(DiffuseSH, DiffuseCoCg, InitialTraceNormal	);
			float[6] SH = IrridianceToSH(DiffuseSHCompute, R);
			TotalSH += vec4(SH[0], SH[1], SH[2], SH[3]);
			TotalCoCg += vec2(SH[4], SH[5]); 
			total_hits ++;
			continue;
		} 
		#endif

		vec3 Normal;
		float Blocktype;

		float T = VoxelTraversalDF(SampledWorldPosition.xyz, R, Normal, Blocktype, false);
		vec3 HitPosition = SampledWorldPosition.xyz + (R * T);

		vec2 UV; 
		vec3 Tangent, Bitangent;
		CalculateVectors(HitPosition, Normal, Tangent, Bitangent, UV); UV.y = 1.0f - UV.y;

		if (T > 0.0f)
		{
			bool ReprojectionSuccessful = false;
			vec2 ScreenSpaceReprojected = vec2(-1.0f);
			
			// Reproject to screen space !

			vec3 Ambient = BaseIndirectDiffuse;
			bool UseLPVGI = false;

			if (u_ReprojectToScreenSpace) {
				vec4 ReprojectedWorldPos;
				ScreenSpaceReprojected = ReprojectReflectionToScreenSpace(HitPosition, Normal, ReprojectionSuccessful, ReprojectedWorldPos);

				if (ReprojectionSuccessful) {
					vec4 ReprojectedSH = texture(u_DiffuseSH, ScreenSpaceReprojected);
					vec2 ReprojectedCoCg = texture(u_DiffuseCoCg, ScreenSpaceReprojected).xy;
					Ambient = SHToIrradianceA(ReprojectedSH, ReprojectedCoCg);
					float ReprojectedVXAO = pow(texture(u_IndirectAO, ScreenSpaceReprojected.xy).x,0.75f);
					if (ReprojectedWorldPos.w>0.0f){
						if (distance(ReprojectedWorldPos.xyz,u_InverseView[3].xyz) < 40) {
							Ambient *= ReprojectedVXAO;
						}
					}
				} 


			}

			if (u_LPVGI && !ReprojectionSuccessful) {
				Ambient = SampleLPV(HitPosition+Normal*0.5f, BaseIndirectDiffuse);
			}


			MaxHitDistance = max(MaxHitDistance, T); Hit = true;
			int reference_id = clamp(int(floor(Blocktype * 255.0f)), 0, 127);

			vec4 texture_ids = vec4(
				float(BlockAlbedoData[reference_id]),
				float(BlockNormalData[reference_id]),
				float(BlockPBRData[reference_id]),
				float(BlockEmissiveData[reference_id])
			);
			
			// I hate this.
			if (reference_id == u_GrassBlockProps[0])
			{
			    if (Normal == NORMAL_LEFT || Normal == NORMAL_RIGHT || Normal == NORMAL_FRONT || Normal == NORMAL_BACK)
				{
					texture_ids.x = u_GrassBlockProps[4];
					texture_ids.y = u_GrassBlockProps[5];
					texture_ids.z = u_GrassBlockProps[6];
				}

				else if (Normal == NORMAL_TOP)
				{
					texture_ids.x = u_GrassBlockProps[1];
					texture_ids.y = u_GrassBlockProps[2];
					texture_ids.z = u_GrassBlockProps[3];
				}

				else if (Normal == NORMAL_BOTTOM)
				{
					texture_ids.x = u_GrassBlockProps[7];
					texture_ids.y = u_GrassBlockProps[8];
					texture_ids.z = u_GrassBlockProps[9];
				}
			}

			mat3 TBN;
			TBN = mat3((Tangent), (Bitangent), (Normal));

			vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(UV,texture_ids.x), 2).rgb;
			vec3 Radiance = SAMPLED_COLOR_MIXED * 0.6f; 
				
			vec4 SampledPBR = textureLod(u_BlockPBRTextures, vec3(UV, texture_ids.z), 3).rgba;
			float AO = pow(SampledPBR.w, 2.0f);

			bool PlayerInShadow = GetPlayerIntersect(HitPosition + Normal*0.035f, NormalizedStrongerDir);
			
			// Compute shadow rays for only 1/4 the reflection samples because performance :p
			if ((ShadowItr < max(SPP / 4, 1))) 
			{
				if (!PlayerInShadow) {
			
					// Try to reuse screen space shadow info to avoid casting a shadow ray : 
					if (ReprojectionSuccessful && u_ReprojectToScreenSpace && InThresholdedScreenSpace(ScreenSpaceReprojected))
					{
						ComputedShadow = texture(u_ShadowTrace, ScreenSpaceReprojected).x;
					}

					else 
					{
						ComputedShadow = GetShadowAt(HitPosition + Normal*0.055f, NormalizedStrongerDir);
					}
				}
				
				else {
				
					ComputedShadow = 1.0f;
				}

				ShadowItr = ShadowItr + 1;
			}

			
			const float AmbientBias = 1.0f;
			Ambient = (Ambient * AmbientBias * clamp(AO, 0.1f, 1.0f)) * vec3(Albedo);

			// Basic normal mapping ->
			// Sampled at a lower LOD
			vec3 NormalMapped = TBN * (textureLod(u_BlockNormalTextures, vec3(UV,texture_ids.y), 3).rgb * 2.0f - 1.0f); 
			
			// Calculate cook torrance brdf ->
			// This is not a 100% accurate, mostly for performance sake (and to reduce variance)
			
			vec3 DirectLighting =  Ambient + CalculateDirectionalLight(HitPosition, 
								   NormalizedStrongerDir, 
								   Radiance, 
								   Albedo, 
								   NormalMapped, 
								   SampledPBR.xyz,
								   ComputedShadow);
			
			if (texture_ids.w > -0.5f) // If the block is not emissive, the read data will be -1!
			{
				float Emissivity = texture(u_BlockEmissiveTextures, vec3(UV, texture_ids.w)).r;
				
				if (Emissivity > 0.2f)
				{
					float m = 20.0f;
					
					// Fix light leak at the edges ->
					float lbiasx = 0.02501f;
                    float lbiasy = 0.03001f;
                    Emissivity *= float(UV.x > lbiasx && UV.x < 1.0f - lbiasx &&
                                 UV.y > lbiasy && UV.y < 1.0f - lbiasy);
					
					// Compute direct lighting 
					DirectLighting = Albedo * max(Emissivity * m, 2.0f);
					
					// Set emissivity mask 
					// This is used as in input to the denoiser
					EmissivityMask = 1.0f;
				}
			}
			
			vec3 Computed;
			Computed = DirectLighting;

			if (u_ReflectPlayer) {
				ComputePlayerReflection(refpos, R, Computed, T);
			}
			
			// Project to 2nd level spherical harmonic 
			float[6] SH = IrridianceToSH(Computed, R);
			TotalSH += vec4(SH[0], SH[1], SH[2], SH[3]);
			TotalCoCg += vec2(SH[4], SH[5]);

			// Store hit distance for reprojection and denoiser ->
			AveragedHitDistance += T; TotalMeaningfulHits += 1.0f;
		}

		else
		{
			vec3 AtmosphereColor;
			GetAtmosphere(AtmosphereColor, R);
			if (u_ReflectPlayer) {
				ComputePlayerReflection(refpos.xyz, R, AtmosphereColor, 10000.0f);
			}
			float[6] SH = IrridianceToSH(AtmosphereColor, R);
			TotalSH += vec4(SH[0], SH[1], SH[2], SH[3]);
			TotalCoCg += vec2(SH[4], SH[5]);
		}


		total_hits++;
	}

	AveragedHitDistance /= max(TotalMeaningfulHits, 0.01f);
	TotalSH /= max(float(total_hits), 1.0f);
	TotalCoCg /= max(float(total_hits), 1.0f);
	
	// Store ->
	o_SH = TotalSH;
	o_CoCg = TotalCoCg;
	o_HitDistance = TotalMeaningfulHits > 0.01f ? AveragedHitDistance : -1.0f;
	o_EmissivityHitMask = EmissivityMask;
	
	// Clamp ->
	o_SH = clamp(o_SH, -100.0f, 100.0f);
	o_CoCg = clamp(o_CoCg, -50.0f, 50.0f);
	o_HitDistance = clamp(o_HitDistance, -10.0f, 200.0f);
	o_EmissivityHitMask = clamp(o_EmissivityHitMask, 0.0f, 1.0f);
}

// Intersection ->

bool IsInVolume(in vec3 pos)
{
    if (pos.x < 0.0f || pos.y < 0.0f || pos.z < 0.0f || 
        pos.x > float(WORLD_SIZE_X - 1) || pos.y > float(WORLD_SIZE_Y - 1) || pos.z > float(WORLD_SIZE_Z - 1))
    {
        return false;    
    }   

    return true;
}

float GetVoxel(ivec3 loc)
{
    if (IsInVolume(loc))
    {
        return texelFetch(u_VoxelData, loc, 0).r;
    }
    
    return 0.0f;
}

float ToConservativeEuclidean(float Manhattan)
{
	return Manhattan == 1 ? 1 : Manhattan * 0.57735026918f;
}

float GetDistance(ivec3 loc)
{
    if (IsInVolume(loc))
    {
         return (texelFetch(u_DistanceFieldTexture, loc, 0).r);
    }
    
    return -1.0f;
}

bool VoxelExists(in vec3 loc)
{
    if (GetVoxel(ivec3(loc)) > 0.0f) 
    {
        return true;
    }

    return false;
}

float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType, bool shadow) 
{
	vec3 initial_origin = origin;
	const float epsilon = 0.01f;
	bool Intersection = false;

	int MinIdx = 0;
	ivec3 RaySign = ivec3(sign(direction));

	int itr = 0;
	int sz = shadow ? 150 : TRACE_LENGTH;

	for (itr = 0 ; itr < sz ; itr++)
	{
		ivec3 Loc = ivec3(floor(origin));
		
		if (!IsInVolume(Loc))
		{
			Intersection = false;
			break;
		}

		float Dist = GetDistance(Loc) * 255.0f; 
		int Euclidean = int(floor(ToConservativeEuclidean(Dist)));

		if (Euclidean == 0)
		{
			break;
		}

		if (Euclidean == 1)
		{
			// Do the DDA algorithm for one voxel 

			ivec3 GridCoords = ivec3(origin);
			vec3 WithinVoxelCoords = origin - GridCoords;
			vec3 DistanceFactor = (((1 + RaySign) >> 1) - WithinVoxelCoords) * (1.0f / direction);

			MinIdx = DistanceFactor.x < DistanceFactor.y && RaySign.x != 0
				? (DistanceFactor.x < DistanceFactor.z || RaySign.z == 0 ? 0 : 2)
				: (DistanceFactor.y < DistanceFactor.z || RaySign.z == 0 ? 1 : 2);

			GridCoords[MinIdx] += RaySign[MinIdx];
			WithinVoxelCoords += direction * DistanceFactor[MinIdx];
			WithinVoxelCoords[MinIdx] = 1 - ((1 + RaySign) >> 1) [MinIdx]; // Bit shifts (on ints) to avoid division

			origin = GridCoords + WithinVoxelCoords;
			origin[MinIdx] += RaySign[MinIdx] * 0.0001f;

			Intersection = true;
		}

		else 
		{
			origin += int(Euclidean - 1) * direction;
		}
	}

	if (Intersection)
	{
		normal = vec3(0.0f);
		normal[MinIdx] = -RaySign[MinIdx];
		blockType = GetVoxel(ivec3(floor(origin)));
		return blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
	}

	return -1.0f;
}


// Ray capsule intersection test ->

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

vec3 capNormal(in vec3 pos, in vec3 a, in vec3 b, in float r)
{
    vec3  ba = b - a;
    vec3  pa = pos - a;
    float h = clamp(dot(pa,ba)/dot(ba,ba),0.0,1.0);
    return (pa - h*ba)/r;
}

bool GetPlayerIntersect(in vec3 WorldPos, in vec3 d)
{
    float x = 0.4;
	vec3 VP = u_ViewerPosition + vec3(-x, -x, +x);
    float t = capIntersect(WorldPos, d, VP, VP + vec3(0.0f, 1.0f, 0.0f), 0.5f);
    return t > 0.0f;
}

void ComputePlayerReflection(in vec3 ro, in vec3 rd, inout vec3 col, float block_t)
{
	float t = capIntersect(ro, rd, u_ViewerPosition - vec3(0.0f, 0.75f, 0.0f), u_ViewerPosition + vec3(0.0f, 0.75f, 0.0f), 0.5f);
	
	if (t > 0.0f)
	{
		if (t < block_t + 0.001f)
		{
			vec3 p = ro + (t * rd);
			vec3 n = capNormal(p, u_ViewerPosition - vec3(0.0f, 0.75f, 0.0f), u_ViewerPosition + vec3(0.0f, 0.75f, 0.0f), 0.5f);
			vec3 albedo = vec3(0.6f);
			float diff = max(dot(n, (u_StrongerLightDirection)), 0.0f);
			col = vec3(diff) * albedo;
			col += vec3(0.150f) * albedo;
		}
	}
}

float GetShadowAt(in vec3 pos, in vec3 ldir)
{
	float T = -1.0f;
	 
	vec3 norm;
	float block;

	if (GetPlayerIntersect(pos, ldir)) { return 1.0f; }
	T = VoxelTraversalDF(pos.rgb, ldir, norm, block, true);

	if (T > 0.0f) 
	{ 
		return 1.0f; 
	}
	
	else 
	{ 
		return 0.0f;
	}
}

float SRGBToLinear(float x){
    return x > 0.04045 ? pow(x * (1 / 1.055) + 0.0521327, 2.4) : x / 12.92;
}

vec3 SRGBToLinearVec3(vec3 x){
    return vec3(SRGBToLinear(x.x),
                SRGBToLinear(x.y),
                SRGBToLinear(x.z));
}

vec3 TemperatureToRGB(float temperatureInKelvins)
{
	vec3 retColor;
	
    temperatureInKelvins = clamp(temperatureInKelvins, 1000, 50000) / 100;
    
    if (temperatureInKelvins <= 66){
        retColor.r = 1;
        retColor.g = clamp01(0.39008157876901960784 * log(temperatureInKelvins) - 0.63184144378862745098);
    } else {
    	float t = temperatureInKelvins - 60;
        retColor.r = clamp01(1.29293618606274509804 * pow(t, -0.1332047592));
        retColor.g = clamp01(1.12989086089529411765 * pow(t, -0.0755148492));
    }
    
    if (temperatureInKelvins >= 66)
        retColor.b = 1;
    else if(temperatureInKelvins <= 19)
        retColor.b = 0;
    else
        retColor.b = clamp01(0.54320678911019607843 * log(temperatureInKelvins - 10) - 1.19625408914);

    return SRGBToLinearVec3(retColor);
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

void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv)
{
	const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
	const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
	const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
	const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
	const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
	const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

    if (CompareVec3(normal, NORMAL_TOP))
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (CompareVec3(normal, NORMAL_BOTTOM))
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (CompareVec3(normal, NORMAL_RIGHT))
    {
        uv = vec2(fract(world_pos.zy));
    }

    else if (CompareVec3(normal, NORMAL_LEFT))
    {
        uv = vec2(fract(world_pos.zy));
    }
    
    else if (CompareVec3(normal, NORMAL_FRONT))
    {
        uv = vec2(fract(world_pos.xy));
    }

     else if (CompareVec3(normal, NORMAL_BACK))
    {
        uv = vec2(fract(world_pos.xy));
    }
}


vec2 ProjectDirection(vec3 Direction, vec2 TextureSize) 
{
	float TileSize = min(floor(TextureSize.x * 0.5f) / 1.5f, floor(TextureSize.y * 0.5f));
	float TileSizeDivided = (0.5f * TileSize) - 1.5f;
	vec2 CurrentCoordinate;

	if (abs(Direction.x) > abs(Direction.y) && abs(Direction.x) > abs(Direction.z)) 
    {
		Direction /= max(abs(Direction.x), 0.001f);
		CurrentCoordinate.x = Direction.y * TileSizeDivided + TileSize * 0.5f;
		CurrentCoordinate.y = Direction.z * TileSizeDivided + TileSize * (Direction.x < 0.0f ? 0.5f : 1.5f);
	} 
    
    else if (abs(Direction.y) > abs(Direction.x) && abs(Direction.y) > abs(Direction.z))
    {
		Direction /= max(abs(Direction.y), 0.001f);
		CurrentCoordinate.x = Direction.x * TileSizeDivided + TileSize * 1.5f;
		CurrentCoordinate.y = Direction.z * TileSizeDivided + TileSize * (Direction.y < 0.0f ? 0.5f : 1.5f);
	} 
    
    else 
    {
		Direction /= max(abs(Direction.z), 0.001f);
		CurrentCoordinate.x = Direction.x * TileSizeDivided + TileSize * 2.5f;
		CurrentCoordinate.y = Direction.y * TileSizeDivided + TileSize * (Direction.z < 0.0f ? 0.5f : 1.5f);
	}

	return CurrentCoordinate / max(TextureSize, 0.01f);
}


vec3 SampleLPVColor(vec3 UV) {
    uint BlockID = texture(u_LPVBlocks, UV).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);
}   


vec3 InterpolateLPVColorDithered(vec3 UV) 
{ 
    const vec3 VolumeResolution = vec3(384.0f, 128.0f, 384.0f);
    vec3 FractTexel = fract(UV * VolumeResolution);
    vec3 LinearOffset = (FractTexel * (FractTexel - 1.0f) + 0.5f) / VolumeResolution;
    vec3 W0 = UV - LinearOffset;
    vec3 W1 = UV + LinearOffset;
    const float DitherWeights[4] = float[4](1.0f, 1.0, 1.0f, 1.0f);
    const float GlobalDitherNoiseWeight = 2.0f;

    vec3 Interpolated = SampleLPVColor(vec3(W0.x, W0.y, W0.z) + LPVDither * DitherWeights[0] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W0.y, W0.z) - LPVDither * DitherWeights[0] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W1.y, W0.z) + LPVDither * DitherWeights[1] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W0.x, W1.y, W0.z) - LPVDither * DitherWeights[1] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W0.x, W1.y, W1.z) + LPVDither * DitherWeights[2] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W1.y, W1.z) - LPVDither * DitherWeights[2] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W0.y, W1.z) + LPVDither * DitherWeights[3] * GlobalDitherNoiseWeight)
		   + SampleLPVColor(vec3(W0.x, W0.y, W1.z) - LPVDither * DitherWeights[3] * GlobalDitherNoiseWeight);
	return max((Interpolated / 8.0), 0.00000001f);
}

vec3 SampleLPVColorTexel(ivec3 Texel, int L) {
    uint BlockID = texelFetch(u_LPVBlocks, Texel, 0).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

// Ground truth triquadratic interpolation, only here for experimentation mostly
vec3 InterpolateLPVColorData(vec3 uv)
{
    vec3 res = vec3(384.0f, 128.0f, 384.0f);
    vec3 q = fract(uv * res);
    ivec3 t = ivec3(uv * res);
    ivec3 e = ivec3(-1, 0, 1);
    vec3 q0 = (q+1.0)/2.0;
    vec3 q1 = q/2.0;	
    vec3 s000 = SampleLPVColorTexel(t + e.xxx, 0);
    vec3 s001 = SampleLPVColorTexel(t + e.xxy, 0);
    vec3 s002 = SampleLPVColorTexel(t + e.xxz, 0);
    vec3 s012 = SampleLPVColorTexel(t + e.xyz, 0);
    vec3 s011 = SampleLPVColorTexel(t + e.xyy, 0);
    vec3 s010 = SampleLPVColorTexel(t + e.xyx, 0);
    vec3 s020 = SampleLPVColorTexel(t + e.xzx, 0);
    vec3 s021 = SampleLPVColorTexel(t + e.xzy, 0);
    vec3 s022 = SampleLPVColorTexel(t + e.xzz, 0);
    vec3 y00 = mix(mix(s000, s001, q0.z), mix(s001, s002, q1.z), q.z);
    vec3 y01 = mix(mix(s010, s011, q0.z), mix(s011, s012, q1.z), q.z);
    vec3 y02 = mix(mix(s020, s021, q0.z), mix(s021, s022, q1.z), q.z);
	vec3 x0 = mix(mix(y00, y01, q0.y), mix(y01, y02, q1.y), q.y);
    vec3 s122 = SampleLPVColorTexel(t + e.yzz, 0);
    vec3 s121 = SampleLPVColorTexel(t + e.yzy, 0);
    vec3 s120 = SampleLPVColorTexel(t + e.yzx, 0);
    vec3 s110 = SampleLPVColorTexel(t + e.yyx, 0);
    vec3 s111 = SampleLPVColorTexel(t + e.yyy, 0);
    vec3 s112 = SampleLPVColorTexel(t + e.yyz, 0);
    vec3 s102 = SampleLPVColorTexel(t + e.yxz, 0);
    vec3 s101 = SampleLPVColorTexel(t + e.yxy, 0);
    vec3 s100 = SampleLPVColorTexel(t + e.yxx, 0);
    vec3 y10 = mix(mix(s100, s101, q0.z), mix(s101, s102, q1.z), q.z);
    vec3 y11 = mix(mix(s110, s111, q0.z), mix(s111, s112, q1.z), q.z);
    vec3 y12 = mix(mix(s120, s121, q0.z), mix(s121, s122, q1.z), q.z);
    vec3 x1 = mix(mix(y10, y11, q0.y), mix(y11, y12, q1.y), q.y);
    vec3 s200 = SampleLPVColorTexel(t + e.zxx, 0);
    vec3 s201 = SampleLPVColorTexel(t + e.zxy, 0);
    vec3 s202 = SampleLPVColorTexel(t + e.zxz, 0);
    vec3 s212 = SampleLPVColorTexel(t + e.zyz, 0);
    vec3 s211 = SampleLPVColorTexel(t + e.zyy, 0);
    vec3 s210 = SampleLPVColorTexel(t + e.zyx, 0);
    vec3 s220 = SampleLPVColorTexel(t + e.zzx, 0);
    vec3 s221 = SampleLPVColorTexel(t + e.zzy, 0);
    vec3 s222 = SampleLPVColorTexel(t + e.zzz, 0);
    vec3 y20 = mix(mix(s200, s201, q0.z), mix(s201, s202, q1.z), q.z);
    vec3 y21 = mix(mix(s210, s211, q0.z), mix(s211, s212, q1.z), q.z);
    vec3 y22 = mix(mix(s220, s221, q0.z), mix(s221, s222, q1.z), q.z);
    vec3 x2 = mix(mix(y20, y21, q0.y), mix(y21, y22, q1.y), q.y);
    return mix(mix(x0, x1, q0.x), mix(x1, x2, q1.x), q.x);
}


vec3 SampleLPVData(vec3 UV)
{    
    UV *= 1.0f/vec3(384.0f,128.0f,384.0f);
	float level = texture(u_LPV, UV).x;
	vec3 InterpolatedColor = vec3(0.0f);

	if (!u_QualityLPVGI) {
		InterpolatedColor = InterpolateLPVColorDithered(UV);
	}

	else {
		InterpolatedColor = InterpolateLPVColorData(UV);
	}

	vec3 FinalInterpolated = vec3(level*100.0f)*pow(pow(InterpolatedColor,vec3(1.0f/2.0f)),vec3(1.0f,1.0f,1.0f/1.40f))*2.0f*vec3(1.,1.,0.9);
	FinalInterpolated = mix(vec3(Luma(FinalInterpolated)), FinalInterpolated, 0.75f);
    return FinalInterpolated;
}

// End 