#version 430 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

layout (location = 0) out float o_Shadow;
layout (location = 1) out float o_IntersectionTransversal; // -> Used as an input to the denoiser 

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler3D u_VoxelData;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2DArray u_AlbedoTextures;
uniform sampler3D u_DistanceFieldTexture;
uniform sampler2D u_BlueNoiseTexture;

uniform bool u_DoFullTrace;
uniform mat4 u_ShadowProjection;
uniform mat4 u_ShadowView;
uniform float u_Time;
uniform float u_FOV;

uniform vec3 u_LightDirection;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform int u_CurrentFrame;

uniform bool u_ContactHardeningShadows;
uniform bool u_ShouldAlphaTest;

uniform vec2 u_Dimensions;
uniform vec2 u_Halton;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
}; 

vec2 g_TexCoords;

vec2 CalculateUV(vec3 world_pos, in vec3 normal);

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

int GetBlockID(float id)
{
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

bool CompareVec3(vec3 v1, vec3 v2) {
	float e = 0.0125f;
	return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
}

bool CompareFloatNormal(float x, float y) {
    return abs(x - y) < 0.02f;
}


float g_K = 0.0f;

bool StopRay(vec3 P, vec3 N, float Type) {
	int BlockID = GetBlockID(Type);
	if (BlockTransparentData[BlockID]==0) {
		return true;
	}

	vec2 SampleUV = CalculateUV(P,N);
	SampleUV.y = 1.0f - SampleUV.y;
	float D = distance(P,u_InverseView[3].xyz);
	int LOD = int(log2(512.0f / (1.0f / D * g_K)));
	float Alpha = textureLod(u_AlbedoTextures, vec3(SampleUV, float(BlockAlbedoData[BlockID])), clamp(LOD-2.0f, 0.0f, 8.0f)).w;
	return Alpha > 0.975f;
}

float VoxelTraversalDF_AlphaTest(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType) 
{
	vec3 initial_origin = origin;
	const float epsilon = 0.01f;
	bool Intersection = false;

	int MinIdx = 0;
	ivec3 RaySign = ivec3(sign(direction));

	int itr = 0;

	for (itr = 0 ; itr < 350 ; itr++)
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
			vec3 tn = vec3(0.0f);
			tn[MinIdx] = -RaySign[MinIdx];
			float bt = GetVoxel(ivec3(floor(origin)));
			if (StopRay(origin, tn, bt)) {
				break;
			}

			else {
				for (int i = 0 ; i < 4 ; i++) {
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
					
					float bt = GetVoxel(ivec3(floor(origin)));
					if (bt > 0.0f) {
						vec3 tn = vec3(0.0f);
						tn[MinIdx] = -RaySign[MinIdx];
						float bt = GetVoxel(ivec3(floor(origin)));
						if (StopRay(origin, tn, bt)) {
							normal = vec3(0.0f);
							normal[MinIdx] = -RaySign[MinIdx];
							blockType = GetVoxel(ivec3(floor(origin)));
							return blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
						}
					}
				}
			}
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

float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType) 
{
	vec3 initial_origin = origin;
	const float epsilon = 0.01f;
	bool Intersection = false;

	int MinIdx = 0;
	ivec3 RaySign = ivec3(sign(direction));

	int itr = 0;

	for (itr = 0 ; itr < 350 ; itr++)
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




vec2 ReprojectShadow (in vec3 pos)
{
	vec4 Projected = u_ShadowProjection * u_ShadowView * vec4(pos, 1.0f);
	Projected.xyz /= Projected.w;
	Projected.xy = Projected.xy * 0.5f + 0.5f;

	return Projected.xy;
}

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


int MIN = -2147483648;
int MAX = 2147483647;
int RNG_SEED;

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


// basic fract(sin) pseudo random number generator
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

const float PI = 3.14159265359f;

const float PHI = 1.61803398874989484820459; 

float gold_noise(in vec2 xy, in float seed)
{
    return fract(tan(distance(xy*PHI, xy)*seed)*xy.x);
}

vec3 SampleCone(vec2 Xi, float CosThetaMax) 
{
    float CosTheta = (1.0 - Xi.x) + Xi.x * CosThetaMax;
    float SinTheta = sqrt(1.0 - CosTheta * CosTheta);
    float phi = Xi.y * PI * 2.0;
    vec3 L;
    L.x = SinTheta * cos(phi);
    L.y = SinTheta * sin(phi);
    L.z = CosTheta;
    return L;
}

//#define USE_GOLD_NOISE
const float GoldenAngle = PI * (3.0 - sqrt(5.0));

vec2 Vogel(uint sampleIndex, uint samplesCount, float Offset)
{
	float r = sqrt(float(sampleIndex) + 0.5f) / sqrt(float(samplesCount));
	float theta = float(sampleIndex) * GoldenAngle + Offset;
	return r * vec2(cos(theta), sin(theta));
}

//#define ANIMATE_BY_GOLDEN_RATIO
//#define USE_EW_WHITE_NOISE


void main()
{
	g_TexCoords = v_TexCoords;
	g_TexCoords += u_Halton * (1.0f / u_Dimensions) * 1.61803398f;

	g_K = 1.0f / (tan(radians(u_FOV) / (2.0f * u_Dimensions.x)) * 2.0f);
	HASH2SEED = (g_TexCoords.x * g_TexCoords.y) * 100.0 * 2.0f;
	HASH2SEED += fract(u_Time) * 102.0f;
	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(800.0f * u_Time);
	hash2(); hash2();

	// Xor shift once!
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
	
	vec4 RayOrigin = GetPositionAt(u_PositionTexture, g_TexCoords).rgba;

	if (RayOrigin.w < 0.0f) {
		o_Shadow = 0.0f;
		o_IntersectionTransversal = 64.;
		return;
	}



	vec3 NormalizedDir = (u_LightDirection);
	vec3 JitteredLightDirection = NormalizedDir;



	if (u_ContactHardeningShadows)
	{
		vec3 RayPosition = RayOrigin.xyz;

		#ifndef USE_EW_WHITE_NOISE 
			#ifdef ANIMATE_BY_GOLDEN_RATIO
				ivec2 TxS = ivec2(textureSize(u_BlueNoiseTexture, 0).xy);
				vec3 Hash = texelFetch(u_BlueNoiseTexture, ivec2(gl_FragCoord.xy*2.0)%TxS, 0).xyz;//texture(u_BlueNoise, g_TexCoords * (u_Dimensions / vec2(TxS / 2.0f))).xy;
				const float GoldenRatio = 1.61803398875f;
				vec3 Xi = mod(Hash + GoldenRatio * (u_CurrentFrame % 240), 1.0f);
			#else
				vec3 Hash;
				int n = u_CurrentFrame % 1024;
				vec2 off = fract(vec2(n * 12664745, n * 9560333) / 16777216.0) * 1024.0;
				ivec2 TextureSize = textureSize(u_BlueNoiseTexture, 0);
				ivec2 SampleTexelLoc = ivec2(gl_FragCoord.xy + ivec2(floor(off))) % TextureSize;
				Hash = texelFetch(u_BlueNoiseTexture, SampleTexelLoc, 0).xyz;
				vec3 Xi = Hash;
			#endif
		#else 
			vec3 Xi = vec3(hash2(), hash2().x);
		#endif
		
		vec3 L = JitteredLightDirection;
		vec3 T = normalize(cross(L, vec3(0.0f, 1.0f, 1.0f)));
		vec3 B = cross(T, L);
		mat3 TBN = mat3(T, B, L);

		// 0.98906604617 -> physically based
		// 0.999825604617
		const float CosTheta = 0.9999505604617f; // -> changed to reduce variance. THIS IS NOT PHYSICALLY CORRECT
		vec3 ConeSample = TBN * SampleCone(Xi.xy, CosTheta); 
        JitteredLightDirection = ConeSample;
	}

	vec3 RayDirection = (JitteredLightDirection);
	vec3 SampledNormal = SampleNormalFromTex(u_NormalTexture, g_TexCoords).rgb;
	
	float NDotL = dot(SampledNormal, RayDirection);

	if (NDotL <= 0.01f) {
		o_Shadow = 1.0f;
		o_IntersectionTransversal = 32.0f / 100.0f; // -> Any transversal. Doesn't matter 
		return;
	}

	vec3 Bias = SampledNormal * vec3(0.06f);
	
	if (true) //(u_DoFullTrace)
	{
		float T = -1.0f;

		float block_at = GetVoxel(ivec3(floor(RayOrigin.xyz + Bias)));
		 
		if (RayOrigin.w > 0.0f) 
		{
			float b; vec3 n;
			T = u_ShouldAlphaTest ? VoxelTraversalDF_AlphaTest(RayOrigin.xyz + Bias, RayDirection, n, b) : VoxelTraversalDF(RayOrigin.xyz + Bias, RayDirection, n, b);
		}

		o_Shadow = float(T > 0.0f || block_at > 0);
		o_IntersectionTransversal = T / 100.0f;
	}
}

vec2 CalculateUV(vec3 world_pos, in vec3 normal)
{
	vec2 uv;

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

	return uv;
}