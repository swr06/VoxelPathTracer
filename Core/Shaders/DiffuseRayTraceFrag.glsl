// Indirect diffuse raytracing 
// We encode the radiance in a 2nd order spherical harmonic 
// To preserve normal map sharpness 
// (It also helps with spatial upscaling!)

#version 450 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

#define USE_COLORED_DIFFUSE // Applies diffuse from the block albedo
#define USE_HEMISPHERICAL_DIFFUSE_SCATTERING 
#define ANIMATE_NOISE // Has to be enabled for temporal filtering to work properly 
#define MAX_VOXEL_DIST 30
#define MAX_BOUNCE_LIMIT 2

// Bayer matrix, used for testing some dithering, unused 
#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

#define EPSILON 0.0001f

// Outputs diffuse indirect
// Specular indirect is handled separately
// Combined based on the fresnel term

layout (location = 0) out vec4 o_SH; // Projected radiance onto the first 2 spherical harmonics.
layout (location = 1) out vec2 o_CoCg; // Stores the radiance color data in YCoCg
layout (location = 2) out float o_Utility; // Using the first SH band to get the luminance is slightly erraneous so I store this. Used as an input for the SVGF denoiser
layout (location = 3) out vec2 o_AOAndSkyLighting; // World Space VXAO (Exaggerate AO at edges, many times better than SSAO) and the amount of skylighting (used for reflection gi + fog)

in vec2 v_TexCoords; // screen space 

in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler3D u_VoxelData;
uniform sampler3D u_DistanceFieldTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_PositionTexture;
uniform samplerCube u_Skymap;

uniform bool CHECKERBOARD_SPP;

uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlockEmissiveTextures;

uniform vec2 u_Dimensions;
uniform vec2 u_Halton;
uniform float u_Time;
uniform bool u_Supersample;

uniform bool u_APPLY_PLAYER_SHADOW;

uniform bool u_UseDirectSampling;

uniform int u_SPP;
uniform int u_CheckerSPP;
uniform int u_CurrentFrame;
uniform int u_CurrentFrameMod512;
uniform int u_CurrentFrameMod128;
uniform bool u_UseBlueNoise;
uniform float u_GISunStrength;
uniform float u_GISkyStrength;

uniform float u_SunVisibility;

uniform vec3 u_ViewerPosition;
uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;

// Temp
uniform mat4 u_ShadowView;
uniform mat4 u_ShadowProjection;
uniform sampler2D u_ShadowMap;
// Temp

uniform float u_DiffuseLightIntensity = 1.0f;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

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
	int sobol_256spp_256d[256 * 256];
	int scramblingTile[128 * 128 * 8];
	int rankingTile[128 * 128 * 8];
};

// Chunk Light Lists 
layout (std430, binding = 4) buffer LightChunkDataSSBO
{
	vec4 LightChunkPositions[];
};

// Offset and size data 
layout (std430, binding = 5) buffer LightChunkDataOffsetsSSBO
{
	ivec2 LightChunkDataOffsets[24*8*24]; // 384 / 16, 128 / 16, 384 / 16
};


vec2 g_TexCoords = vec2(0.0f);


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



// Function prototypes
float nextFloat(inout int seed, in float min, in float max);
float nextFloat(inout int seed, in float max);
vec3 cosWeightedRandomHemisphereDirection(const vec3 n); 
float nextFloat(inout int seed); 
int nextInt(inout int seed); 
vec3 GetSkyColorAt(vec3 rd);
vec3 RandomPointInUnitSphereRejective();
vec3 CosineSampleHemisphere_2(float u1, float u2);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
float GetShadowAt(vec3 pos, in vec3 ldir);
vec2 ReprojectShadow(vec3);
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv);
float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType, in int dist);
bool voxel_traversal(vec3 origin, vec3 direction, inout float block, out vec3 normal, out vec3 world_pos, int dist);
vec3 GetBRDF(vec3 normal, vec3 incident, vec3 sunDir, vec3 albedo, vec2 PBR);
float DiffuseHammon(vec3 normal, vec3 viewDir, vec3 lightDir, float roughness);
vec3 CreateDiffuseRay(in vec3 world_pos, vec3 radiance, in vec3 albedo, in float shadow, float NDotL, vec3 PBR, vec3 N, vec3 rD, vec3 sdir, out vec3 weight);

// Globals
vec3 g_Normal;
int RNG_SEED = 0;

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

const bool CAUSTICS = false;

// Colors 
const vec3 SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * (16.0f);
const vec3 NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 1.5f; 
const vec3 DUSK_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.9f; 

// basic fract(sin) pseudo random number generator
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

float hash1() 
{
	return fract(sin(HASH2SEED += 0.1) * 43758.5453123);
}

// -----

struct MISPDF {
	float SphericalPDFSample; // Approximated PDF using solid angles and distances
	float LightChoosePDF; // 1.0f / C
};

float lineDistance(vec2 a, vec2 b, vec2 p) 
{
    vec2 pa = p-a;
    vec2 ba = b-a;
	float t = clamp(dot(pa,ba)/dot(ba,ba), 0.0, 1.0);
    return length(pa-ba*t);
}

void basis(in vec3 n, out vec3 b1, out vec3 b2)
{
    float sign_ = sign(n.z);
	float a = -1.0 / (sign_ + n.z);
	float b = n.x * n.y * a;
	b1 = vec3(1.0 + sign_ * n.x * n.x * a, sign_ * b, -sign_ * n.x);
	b2 = vec3(b, sign_ + n.y * n.y * a, -n.y);
}

vec3 localToWorld(in vec3 localDir, in vec3 normal) 
{
    vec3 a, b;
    basis( normal, a, b );
	return localDir.x*a + localDir.y*b + localDir.z*normal;
}

float PdfWtoA(float aPdfW, float aDist2, float aCosThere)
{
    if(aDist2 < EPSILON)
	{
        return 0.0;
	}

    return aPdfW * abs(aCosThere) / aDist2;
}

float PdfAtoW(float aPdfA, float aDist2, float aCosThere)
{
    float absCosTheta = abs(aCosThere);

    if(absCosTheta < EPSILON)
	{
        return 0.0;
    }

    return aPdfA * aDist2 / absCosTheta;
}

vec2 uniformPointWithinCircle(in float radius, in float Xi1, in float Xi2) 
{
    float r = radius * sqrt(1.0 - Xi1);
    float theta = Xi2 *2.0f * PI;
	return vec2(r * cos(theta), r * sin(theta));
}

vec3 uniformDirectionWithinCone(in vec3 d, in float phi, in float sina, in float cosa) 
{    
	vec3 w = normalize(d);
    vec3 u = normalize(cross(w.yzx, w));
    vec3 v = cross(w, u);
	return (u * cos(phi) + v * sin(phi)) * sina + w * cosa;
}

float misWeightPower( in float a, in float b ) {
    float a2 = a*a;
    float b2 = b*b;
    float a2b2 = a2 + b2;
    return a2 / a2b2;
}

float misWeightBalance(in float a, in float b) 
{
    float ab = a + b;
    return a / ab;
}

//#define MIS_HEURISTIC_POWER

float misWeight(in float pdfA, in float pdfB) {
#ifdef MIS_HEURISTIC_POWER
    return misWeightPower(pdfA,pdfB);
#else
    return misWeightBalance(pdfA,pdfB);
#endif
}


// Xi -> Sample points
vec3 SampleCone(vec2 Xi, float CosThetaMax) 
{
    float CosTheta = (1.0f - Xi.x) + Xi.x * CosThetaMax;
    float SinTheta = sqrt(1.0f - CosTheta * CosTheta);
    float phi = Xi.y * PI * 2.0f;
    vec3 L;
    L.x = SinTheta * cos(phi);
    L.y = SinTheta * sin(phi);
    L.z = CosTheta;
    return L;
}

float DistanceSquared(vec3 P1, vec3 P2) {
	vec3 X = abs(P1 - P2);
	return dot(X,X);
}

float LambertianPDF(float CosTheta) {
	return CosTheta / PI; // Simplest form, dot(n,r)/pi 
}

// Generate PDF based on solid angle
float SphericalLightPDF(in vec3 x, in vec3 wi, float d, in vec3 n1, vec3 sO, float sR2)
{
    float SolidAngle;
    vec3 w = sO - x; //direction to light center
	float dc_2 = dot(w, w); //squared distance to light center
    float dc = sqrt(dc_2); //distance to light center
    
    if(dc_2 > sR2) 
	{
    	float sin_theta_max_2 = clamp(sR2 / dc_2, 0.0, 1.0);
		float cos_theta_max = sqrt( 1.0 - sin_theta_max_2 );
    	SolidAngle = (2.0f * PI) * (1.0 - cos_theta_max);
    }
	
	else 
	{ 
    	SolidAngle = (4.0f * PI);
    }
    
    return 1.0f / SolidAngle;
}


struct LightSamplingRecord
{
    vec3 w;
    float d;
    float pdf;
};

void SampleLight(vec3 x, vec3 P, vec2 Xi, out LightSamplingRecord sampleRec) 
{
	float Xi1 = Xi.x;
	float Xi2 = Xi.y;
	float sph_r2 = 1.0f;
    vec3 sph_p = P;
    vec3 w = sph_p - x;		
	float dc_2 = dot(w, w);	
    float dc = sqrt(dc_2);	
    float sin_theta_max_2 = sph_r2 / dc_2;
	float cos_theta_max = sqrt( 1.0 - clamp( sin_theta_max_2, 0.0, 1.0 ) );
    float cos_theta = mix(cos_theta_max, 1.0, Xi1);
    float sin_theta_2 = 1.0 - cos_theta*cos_theta;
    float sin_theta = sqrt(sin_theta_2);
    sampleRec.w = uniformDirectionWithinCone( w, (2.0f*PI)*Xi2, sin_theta, cos_theta );
    sampleRec.pdf = 1.0/( (2.0f*PI) * (1.0 - cos_theta_max) );
    sampleRec.d = dc*cos_theta - sqrt(sph_r2 - dc_2*sin_theta_2);
}

float signum(float x) {
	return uintBitsToFloat((floatBitsToUint(x) & 0x80000000u) | floatBitsToUint(1.0));
}

mat3 GenerateRotationMatrix(vec3 from, vec3 to) 
{
	float cosine = dot(from, to);

	float tmp = signum(cosine);
	tmp = 1.0 / (tmp + cosine);

	vec3 axis = cross(to, from);
	vec3 tmpv = axis * tmp;

	return mat3(
		axis.x * tmpv.x + cosine, axis.x * tmpv.y - axis.z, axis.x * tmpv.z + axis.y,
		axis.y * tmpv.x + axis.z, axis.y * tmpv.y + cosine, axis.y * tmpv.z - axis.x,
		axis.z * tmpv.x - axis.y, axis.z * tmpv.y + axis.x, axis.z * tmpv.z + cosine
	);
}

vec2 ConcentricSampleDisk(vec2 random) 
{
    random = 2.0 * random - 1.0;
    if (all(equal(random, vec2(0)))) return vec2(0);

    float theta, r;

    if (abs(random.x) > abs(random.y)) 
	{
        r = random.x;
        theta = (PI / 4.0f) * (random.y / random.x);
    }
	
	else 
	{
        r = random.y;
        theta = (PI / 2.0f) - (PI / 4.0f) * (random.x / random.y);
    }
    return r * vec2(cos(theta), sin(theta));
}

// Stochastically selects a light from the light chunks 
// s : boolean valid, P : sample position

vec3 StochasticallySampleLight(vec3 P, out bool s, out float LightSelectPDF, out ivec3 LightPositionFetch) 
{
	s = true;

	LightSelectPDF = 1.0f;
	ivec3 block_loc = ivec3(floor(P));
	
	// 16 ^ 3 chonks.
	int cx = int(floor(float(block_loc.x) / float(16)));
	int cy = int(floor(float(block_loc.y) / float(16)));
	int cz = int(floor(float(block_loc.z) / float(16)));

	int OffsetArrayFetchLocation = (cz * 24 * 8) + (cy * 24) + cx;

	ivec2 ChunkData = LightChunkDataOffsets[OffsetArrayFetchLocation];

	if (ChunkData.y <= 0 || ChunkData.x <= 0) {
		s = false;
		return vec3(0.0f);
	}

	s = true;
	int Range = abs(ChunkData.y);


	vec3 BestPosition = vec3(1000.0f);
	float BestDistance = 1000.0f;

	float IterationMixer = 1.0f - exp(-(float(Range) / (16.0f*16.0f)));
	int Iterations = clamp(int(floor(mix(3.0f, 24.0f, IterationMixer))),3,24);

	float Hash = hash1();
	
	float SmoothCutoff = mix(2.0f, 7.0f, Hash);

	for (int Sample = 0 ; Sample < Iterations ; Sample++) { 

		vec2 xi = hash2();
		int RandomLight = int(floor(mix(0.0f, float(Range), xi.x)));
		RandomLight = clamp(RandomLight, 0, Range);
		vec3 RandomPosition = LightChunkPositions[ChunkData.x + RandomLight].xyz;

		float d2 = DistanceSquared(P, RandomPosition);

		float Error = abs(d2 - BestDistance);

		if (Error > SmoothCutoff) {

			BestDistance = d2;
			BestPosition = RandomPosition;
		}

		else if (Error < 32.0f && xi.y > 0.97f) // 3% chance of selecting a far ish away light 
		{
			BestDistance = d2;
			BestPosition = RandomPosition;
		}

	}

	LightPositionFetch = ivec3(floor(BestPosition));

	if (abs(distance(BestPosition,P)) > 24.0f || BestDistance > 64.0f) {
		s = false;
		return vec3(0.0f);
	}

	// Todo : Use a CDF to get an accurate probability density function

	// Approximate light selection pdf 
	const float SelectionBias = 5.0f;

	LightSelectPDF = min(Range, SelectionBias) / max(float(Range), 1.0f); 
	LightSelectPDF = clamp(LightSelectPDF, 0.001f, 1.0f);
	
	s = true;

	return BestPosition;
}


vec3 SampleLightDisk(vec3 P, vec3 Center, vec3 Normal, vec2 Xi, out float DiskPDF) {
	
	// Ideal constants for voxels 
	const float DiskRadius = 0.45f; 
	const float Area = PI / 2.0f;

	// Sample concentric disk 
	vec2 DiskSample = ConcentricSampleDisk(Xi);

	vec3 ScaledPoint = vec3(DiskSample.x, 0.0, DiskSample.y) * DiskRadius;
	vec3 DiskDirection = normalize(P - Center);
	mat3 RotationMatrix = GenerateRotationMatrix(vec3(0,1,0), DiskDirection);

	// Rotate Point 
	ScaledPoint = RotationMatrix * ScaledPoint;

	// Translate 
	ScaledPoint += Center;

	vec3 ScaledPointToOriginDir = normalize(ScaledPoint - P);
	float Distance = clamp(distance(P, ScaledPoint), 0.001, 99.0);
	
	// Generate Disk based on surface area :
	DiskPDF = (0.45f + Distance * Distance) / max(Area * dot(DiskDirection, -ScaledPointToOriginDir) * dot(Normal, ScaledPointToOriginDir), 0.001f);
	
	return ScaledPointToOriginDir;
}

//#define SAMPLE_CONE_DL

vec3 GetStochasticDirectLightDir(bool DontGuide, vec3 O, vec3 N, inout float PDF, inout bool SampleSuccess, inout float sPDF, out ivec3 LightPositionFetch) {
	
	// Sample lambert : 
	if (DontGuide) {
		vec3 Reflected =  cosWeightedRandomHemisphereDirection(N);
		PDF = 1.0f;
		sPDF = dot(N, Reflected) / PI;
	
		return Reflected;
	}
	
	bool Success = false;
	vec3 SelectedLight = StochasticallySampleLight(O, Success, PDF, LightPositionFetch);

	SampleSuccess = Success;

	if (!Success) {
		
		vec3 Reflected =  cosWeightedRandomHemisphereDirection(N);
		PDF = 1.0f;
		sPDF = dot(N, Reflected) / PI;

		return Reflected;
	}

	#ifdef SAMPLE_CONE_DL
		const vec3 Basis = vec3(0.0f, 1.0f, 1.0f);
		vec3 L = normalize(SelectedLight - O);
		vec3 T = normalize(cross(L, Basis));
		vec3 B = cross(T, L);
		mat3 TBN = mat3(T, B, L);
		const float CosTheta = 0.985f; 
		vec3 ConeSample = TBN * SampleCone(hash2(), CosTheta);
		return ConeSample;
	#else 
		// Sample a perfect sized disk ->

		vec3 SampleDirection = SampleLightDisk(O + N * 0.07f, SelectedLight, N, hash2(), sPDF);
		return SampleDirection;
	#endif
}

// ---

// Simplified diffuse brdf
vec3 SunBRDF(in vec3 world_pos, vec3 radiance, in vec3 albedo, in float shadow, float NDotL, vec3 PBR, vec3 N, vec3 rD, vec3 sdir)
{
	if (!CAUSTICS) {
		return albedo * DiffuseHammon(N, -rD, sdir, PBR.x) * (radiance*3.5f) * (1.0f - shadow) * PI;
	}

	else {
		return (GetBRDF(N, -rD, sdir, albedo, PBR.xy) * (radiance*4.0f)) * (1.0f - shadow);
	}
} 

vec3 SunBRDFWithoutRadiance(in vec3 world_pos, vec3 radiance, in vec3 albedo, in float shadow, float NDotL, vec3 PBR, vec3 N, vec3 rD, vec3 sdir)
{
	if (!CAUSTICS) {
		return albedo * DiffuseHammon(N, -rD, sdir, PBR.x) * (1.0f - shadow) * PI;
	}

	else {
		return GetBRDF(N, -rD, sdir, albedo, PBR.xy) * (1.0f - shadow);
	}
} 


vec3 LIGHT_COLOR;
vec3 StrongerLightDirection;
const bool g_AverageBounces = false;
float EmissivityMultiplier = 0.0f;

vec3 DiffuseRayBRDF(vec3 N, vec3 I, vec3 D, float Roughness)
{
    vec3 BRDF = vec3(1.0f) * DiffuseHammon(N, -I, D, Roughness);
    return BRDF;
}

bool Moonstronger=false;

vec4 CalculateDiffuse(in vec3 initial_origin, in vec3 input_normal, out vec3 odir, out bool Skyhit)
{
	Skyhit = false;
	
	float bias = 0.06f;
	//float LightSamplePDF = 1.0f; // x/lc (where x is a strength value and lc is the light count)
	//float SamplingPDF = 1.0f;
	//bool SampleSuccessful = false;
	//vec3 Li = GetStochasticDirectLightDir(!GuideSample, initial_origin, input_normal, LightSamplePDF, SampleSuccessful, SamplingPDF);
	Ray new_ray = Ray(initial_origin + input_normal * bias, cosWeightedRandomHemisphereDirection(input_normal));


	float ao = 1.0f;

	vec3 RayContribution = vec3(0.0f);
	vec3 RayThroughput = vec3(1.0f);
	vec3 SummedContribution_t = vec3(0.0f);

	for (int i = 0 ; i < MAX_BOUNCE_LIMIT ; i++)
	{
		if (i == 0) {
			odir = new_ray.Direction;
		}

		vec3 tangent_normal;
		
		vec3 HitNormal; 
		float HitBlock;
		float T = VoxelTraversalDF(new_ray.Origin, new_ray.Direction, HitNormal, HitBlock, MAX_VOXEL_DIST);
		int tex_ref = clamp(int(floor(HitBlock * 255.0f)), 0, 127);
		bool Intersect = T > 0.0f;
		vec3 IntersectionPosition = new_ray.Origin + (new_ray.Direction * T);

		if (Intersect && HitBlock > 0) 
		{ 
			vec2 txc; 
			CalculateUV(IntersectionPosition, HitNormal, txc);
			vec2 TextureIndexes = vec2(
				float(BlockAlbedoData[tex_ref]),
				float(BlockEmissiveData[tex_ref])
			);

			vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(txc, TextureIndexes.r), 3.0f).rgb; // 512, 256, 128, 64, 32, 16
			vec3 PBR = textureLod(u_BlockPBRTextures, vec3(txc, TextureIndexes.r), 2.0f).rgb; // 512, 256, 128, 64, 32, 16
			float Emmisivity = 0.0f;

			if (TextureIndexes.y >= 0.0f)
			{
				float SampledEmmisivity = texture(u_BlockEmissiveTextures, vec3(txc, TextureIndexes.y)).r;
				Emmisivity = SampledEmmisivity * EmissivityMultiplier * u_DiffuseLightIntensity;
			}

			float NDotL = max(dot(HitNormal, StrongerLightDirection), 0.0f);
			vec3 bias_shadow = (HitNormal * 0.045f);
			float ShadowAt;

			// Don't cast shadow rays at night as the moon doesn't contribute to gi much here
			// Cast more samples instead to reduce the variance
			if (Moonstronger)
			{ 
				ShadowAt = 1.0f;
			}

			else {
				ShadowAt = (NDotL < 0.001f) ? 0.0f : GetShadowAt(IntersectionPosition + bias_shadow, StrongerLightDirection);
			}

			//radiance += (Emmisivity * Albedo)+BRDF;
			
			vec3 EmmisivityColor = (Emmisivity * mix(1.0,1.0,float(u_SunVisibility)) * Albedo);
			vec3 SUNBRDF = SunBRDF(IntersectionPosition, LIGHT_COLOR, Albedo, ShadowAt, NDotL, PBR, HitNormal, new_ray.Direction, StrongerLightDirection);

			vec3 NewDirection = vec3(0.0f);
			NewDirection = cosWeightedRandomHemisphereDirection(HitNormal); // create a direction based on the lambertian diffuse brdf
			float CosTheta = clamp(dot(HitNormal, NewDirection), 0.0f, 1.0f); 
			float PDF = max(CosTheta / PI, 0.00001f); // lambert brdf pdf -> dot(n,r)/pi 
			
			// Weight by Lambertian diffuse brdf because i use a cosine weighted hemisphere function
			vec3 Attenuation = DiffuseRayBRDF(HitNormal, new_ray.Direction, NewDirection, PBR.x); 

			RayContribution += RayThroughput * SUNBRDF; // Sun GI
			RayContribution += EmmisivityColor * RayThroughput; // Emissive GI
			RayThroughput *= Attenuation / PDF; 

			SummedContribution_t += EmmisivityColor + SUNBRDF;

			new_ray.Direction = NewDirection;
			new_ray.Origin = IntersectionPosition + HitNormal * bias;
		} 

		else 
		{	
			float x = mix(1.0, 1.05f, u_SunVisibility);
			x = clamp(x*1.0f*u_GISkyStrength,0.0f,5.0f);
			vec3 sky =  (GetSkyColorAt(new_ray.Direction) * x);
			RayContribution += sky * RayThroughput;
			SummedContribution_t += sky;
			Skyhit = true;
			break;
		}


		if (i == 0)
		{
			const bool band_ao = false;
			const float dao = band_ao ? 3.5f : 2.0f;

			if (!band_ao) {

				if (T < dao && T > 0.0f) 
				{
					// Calculate ao on first bounce
					
					ao = max(T / dao, 0.0f);
				}
			}

			else {

				if (T < dao && T > 0.0f) 
				{
					ao = 1.0f - float(T*T<0.225f);
				}
			}
		}

	}

	return g_AverageBounces ? vec4(SummedContribution_t / float(MAX_BOUNCE_LIMIT), ao) : vec4(RayContribution, ao); 
}

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

vec3 DirectSample(in vec3 initial_origin, in vec3 input_normal, out vec3 odir)
{
	float bias = 0.06f;

	// PDF of the stochastic light sample ->
	float LightSamplePDF = 1.0f;

	// PDF of the disk sample ->
	float SamplingPDF = 1.0f;
	
	// Success?
	bool SampleSuccessful = false;

	// The light chosen to sample 
	ivec3 LightPosition = ivec3(0);

	// Normalized sample direction ->
	vec3 Wi = GetStochasticDirectLightDir(false, initial_origin, input_normal, LightSamplePDF, SampleSuccessful, SamplingPDF, LightPosition);
	odir = Wi;

	if (!SampleSuccessful) {
		return vec3(0.0f);
	}

	Ray new_ray = Ray(initial_origin + input_normal * bias, Wi);


	vec3 HitNormal; 
	float HitBlock;

	// Intersect Shadow ray ->
	float T = VoxelTraversalDF(new_ray.Origin, new_ray.Direction, HitNormal, HitBlock, MAX_VOXEL_DIST);

	// Intersection 
	int tex_ref = clamp(int(floor(HitBlock * 255.0f)), 0, 127); 
	bool Intersect = T > 0.0f;
	vec3 IntersectionPosition = new_ray.Origin + (new_ray.Direction * T);


	vec3 Lo = vec3(0.0f);

	if (Intersect && HitBlock > 0) 
	{ 
		vec2 txc; 
		CalculateUV(IntersectionPosition, HitNormal, txc);
		vec2 TextureIndexes = vec2(
			float(BlockAlbedoData[tex_ref]),
			float(BlockEmissiveData[tex_ref])
		);

		vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(txc, TextureIndexes.r), 3.0f).rgb; // 512, 256, 128, 64, 32, 16
		float Emmisivity = 0.0f;

		if (TextureIndexes.y >= 0.0f)
		{
			float SampledEmmisivity = texture(u_BlockEmissiveTextures, vec3(txc, TextureIndexes.y)).r;
			Emmisivity = SampledEmmisivity * EmissivityMultiplier * u_DiffuseLightIntensity;
		}

		// Fetch incoming radiance ->
		vec3 Li = (Emmisivity * mix(1.0,1.0,float(u_SunVisibility)) * Albedo);
		
		float CosTheta = clamp(dot(HitNormal, Wi), 0.0f, 1.0f); 
		float LambertianPDF = max(CosTheta / PI, 0.00001f);
		const vec3 LambertBRDF = vec3(1.0f / PI);  // Lambert brdf 

		float EvaluatedWeight = 1.0f;

		float DiskPDF = SamplingPDF;
		float DiskPDFSquared = DiskPDF * DiskPDF;

		//EvaluatedWeight = misWeight(DiskPDF, LambertianPDF);

		Lo += Li * (1.0f / DiskPDF) * PI; //(Li * weightDisk * (1.0f / DiskPDF));
	} 


	return Lo; 
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

// from quake2rtx project
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

bool CompareFloatNormal(float x, float y) {
    return abs(x - y) < 0.02f;
}

vec3 GetNormalFromID(float n) {
	const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f));
    int idx = int(round(n*10.0f));

    if (idx > 5) {
        return vec3(0.5f);
    }

    return Normals[idx];
}

vec3 SampleNormal(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}

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
	return Noise;
}

void main()
{
	// Compute globals

	g_TexCoords = v_TexCoords;

	if (u_Supersample) {
		g_TexCoords += (u_Halton * 0.75f) / u_Dimensions;
	}

	bool SunStronger = -u_SunDirection.y < 0.01f ? true : false;
	float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 2.9), 0.0f, 1.0f);
    vec3 SunColor = mix(SUN_COLOR, DUSK_COLOR, DuskVisibility);
	LIGHT_COLOR = SunStronger ? SunColor : NIGHT_COLOR;
	LIGHT_COLOR *= 0.4f*u_GISunStrength;
	StrongerLightDirection = SunStronger ? u_SunDirection : u_MoonDirection;
	Moonstronger = !SunStronger;
	EmissivityMultiplier = Moonstronger ? 15.0f : 12.0f;
	// 


	o_AOAndSkyLighting = vec2(1.0f, 0.0f);
	

    #ifdef ANIMATE_NOISE
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(fract(u_Time) * 1000);
	#else
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x);
	#endif

	// XOR shift
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;

	// Animate RNG
	HASH2SEED = (g_TexCoords.x * g_TexCoords.y) * 64.0 * 8.0f;
	HASH2SEED += fract(u_Time) * 128.0f;

	vec4 Position = GetPositionAt(u_PositionTexture, g_TexCoords);
	vec3 Normal = GetNormalFromID(texture(u_NormalTexture, g_TexCoords).x);

	o_Utility = GetLuminance(vec3(0.0f)); 


	if (Position.w < 0.0f)
	{
		float SH[6] = IrridianceToSH(texture(u_Skymap, normalize(v_RayDirection)).xyz * 2.66f, Normal);
		o_SH = vec4(SH[0], SH[1], SH[2], SH[3]);
		o_CoCg.xy = vec2(SH[4], SH[5]);
		return;
	}

	float AccumulatedAO = 0.0f;

	int SPP = clamp(u_SPP, 1, 32);

	if (CHECKERBOARD_SPP) {
		bool CheckerStep = int(gl_FragCoord.x + gl_FragCoord.y) % 2 == u_CurrentFrame % 2;
		SPP = int(mix(u_SPP, u_CheckerSPP, float(CheckerStep)));
	}

	SPP = clamp(SPP, 1, 32);

	vec4 TotalSHy = vec4(0.0f);
	vec2 CoCg = vec2(0.0f);
	vec3 radiance = vec3(0.0f);
	float Skyhits = 0.0f;

	if (Moonstronger) {
		SPP*=2; // We DONT cast shadow rays if the moon is stronger, so we save two shadow rays per diffuse ray, double the spp to make up for this
	}
	
	
	for (int s = 0 ; s < SPP ; s++)
	{
		vec3 DirectSampleDir;
		bool ss=false;
	    //vec3 DirectSample = DirectSample(Position.xyz, Normal, DirectSampleDir);
		vec3 d = vec3(0.0f);
		vec4 x = u_UseDirectSampling ? vec4(DirectSample(Position.xyz, Normal, DirectSampleDir), 1.0f) : CalculateDiffuse(Position.xyz, Normal, d, ss);
		x.xyz = clamp(x.xyz,0.0f,8.0f);

		radiance += x.xyz;
		
		AccumulatedAO += x.w;
		
		float SH[6] = IrridianceToSH(x.xyz, d);
		TotalSHy += vec4(SH[0], SH[1], SH[2], SH[3]);
		CoCg += vec2(SH[4], SH[5]);

		Skyhits += float(ss);
	}


	AccumulatedAO /= SPP;
	TotalSHy /= SPP;
	CoCg /= SPP;
	radiance /= SPP;
	Skyhits /= SPP;


	o_SH = TotalSHy;
	o_CoCg.xy = CoCg;
	o_Utility = max(GetLuminance(radiance),0.01f);
	o_AOAndSkyLighting.x = AccumulatedAO;
	o_AOAndSkyLighting.y = Skyhits;

	// Clamp everything
	o_AOAndSkyLighting = clamp(o_AOAndSkyLighting,0.0f,1.0f);

	o_SH = clamp(o_SH, -100.0f, 100.0f);
	o_CoCg = clamp(o_CoCg, -100.0f, 100.0f);
	o_Utility = clamp(o_Utility, 0.001f, 64.0f);
}

vec3 lerp(vec3 v1, vec3 v2, float t)
{
	return (1.0f - t) * v1 + t * v2;
}


// RNG
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

// RNG

vec3 cosWeightedRandomHemisphereDirection(const vec3 n) 
{
  	vec2 r = vec2(0.0f);

	if (!u_UseBlueNoise) {
		r = vec2(hash2());
	} 
	
	else {
		r = SampleBlueNoise2D(u_CurrentFrameMod128);
	}

	float PI2 = 2.0f * PI;
	vec3  uu = normalize(cross(n, vec3(0.0,1.0,1.0)));
	vec3  vv = cross(uu, n);
	float ra = sqrt(r.y);
	float rx = ra * cos(PI2 * r.x); 
	float ry = ra * sin(PI2 * r.x);
	float rz = sqrt(1.0 - r.y);
	vec3  rr = vec3(rx * uu + ry * vv + rz * n );
    
    return normalize(rr);
}

// vec3 SunBRDF(in vec3 world_pos, vec3 radiance, in vec3 albedo, in float shadow, float NDotL, vec3 PBR, vec3 N, vec3 rD, vec3 sdir)
vec3 GGX_LIGHT(vec3 n, vec3 v, vec3 l, vec2 RM, vec3 albedo) ;

vec3 CreateDiffuseRay(in vec3 world_pos, vec3 radiance, in vec3 albedo, in float shadow, float NDotL, vec3 PBR, vec3 N, vec3 rD, vec3 sdir, out vec3 weight)
{
	const float OneOverPI = 1.0f / PI;
    vec3 DiffuseDirection;
	float PDF;
    DiffuseDirection = cosWeightedRandomHemisphereDirection(N);
    float CosTheta = clamp(dot(N, DiffuseDirection), 0.0f, 1.0f);
    PDF = max(CosTheta * OneOverPI, 0.00001f); // dot(n,r)/pi
    weight = DiffuseRayBRDF(N, rD, DiffuseDirection, PBR.x) / PDF;
    return DiffuseDirection;
}

vec3 GetSkyColorAt(vec3 rd) 
{
	rd.y = clamp(rd.y, 0.125f, 1.5f);
    return texture(u_Skymap, (rd)).rgb;
}


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

float GetManhattanDist(vec3 p1, vec3 p2)
{
	float Manhattan = abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
	return Manhattan;
}

float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType, in int dist) 
{
	vec3 initial_origin = origin;
	const float epsilon = 0.01f;
	bool Intersection = false;

	int MinIdx = 0;
	ivec3 RaySign = ivec3(sign(direction));

	int itr = 0;

	for (itr = 0 ; itr < dist ; itr++)
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

bool voxel_traversal(vec3 origin, vec3 direction, inout float block, out vec3 normal, out vec3 world_pos, int dist)
{
	const vec3 BLOCK_CALCULATED_NORMALS[6] = vec3[]
	(
			vec3(1.0, 0.0, 0.0),
			vec3(-1.0, 0.0, 0.0),
			vec3(0.0, 1.0, 0.0),
			vec3(0.0, -1.0, 0.0),
			vec3(0.0, 0.0, 1.0),
			vec3(0.0, 0.0, -1.0)
	);
	
	world_pos = origin;

	vec3 Temp;
	vec3 VoxelCoord; 
	vec3 FractPosition;

	Temp.x = direction.x > 0.0 ? 1.0 : 0.0;
	Temp.y = direction.y > 0.0 ? 1.0 : 0.0;
	Temp.z = direction.z > 0.0 ? 1.0 : 0.0;

	vec3 plane = floor(world_pos + Temp);

	for (int x = 0; x < dist; x++)
	{
		if (!IsInVolume(world_pos))
		{
			break;
		}

		vec3 Next = (plane - world_pos) / direction;
		int side = 0;

		if (Next.x < min(Next.y, Next.z)) 
		{
			world_pos += direction * Next.x;
			world_pos.x = plane.x;
			plane.x += sign(direction.x);
			side = 0;
		}

		else if (Next.y < Next.z) 
		{
			world_pos += direction * Next.y;
			world_pos.y = plane.y;
			plane.y += sign(direction.y);
			side = 1;
		}

		else 
		{
			world_pos += direction * Next.z;
			world_pos.z = plane.z;
			plane.z += sign(direction.z);
			side = 2;
		}

		VoxelCoord = (plane - Temp);

		int Side = ((side + 1) * 2) - 1;

		if (side == 0) 
		{
			if (world_pos.x - VoxelCoord.x > 0.5)
			{
				Side = 0;
			}
		}

		else if (side == 1)
		{
			if (world_pos.y - VoxelCoord.y > 0.5)
			{
				Side = 2;
			}
		}

		else 
		{
			if (world_pos.z - VoxelCoord.z > 0.5)
			{
				Side = 4;
			}
		}

		normal = BLOCK_CALCULATED_NORMALS[Side];
		block = GetVoxel(ivec3(VoxelCoord.xyz));

		if (block > 0)
		{
			return true; 
		}
	}

	return false;
}


/*

// initial DDA algorithm 
// thanks to telo for the help with understanding it 


// Project to cube to make sure that the ray always starts at the volume -> 
float ProjectToCube(vec3 ro, vec3 rd) 
{	
	float tx1 = (0 - ro.x) / rd.x;
	float tx2 = (MapSize.x - ro.x) / rd.x;

	float ty1 = (0 - ro.y) / rd.y;
	float ty2 = (MapSize.y - ro.y) / rd.y;

	float tz1 = (0 - ro.z) / rd.z;
	float tz2 = (MapSize.z - ro.z) / rd.z;

	float tx = max(min(tx1, tx2), 0);
	float ty = max(min(ty1, ty2), 0);
	float tz = max(min(tz1, tz2), 0);

	float t = max(tx, max(ty, tz));
	
	return t;
}

// My *FIRST* ever DDA implementation, amazing. I know.
float DDA(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType) 
{
	vec3 origin = orig;
	const float epsilon = 0.001f;
	float t1 = max(ProjectToCube(origin, direction) - epsilon, 0.0f);
	origin += t1 * direction;

	int mapX = int(floor(origin.x));
	int mapY = int(floor(origin.y));
	int mapZ = int(floor(origin.z));

	float sideDistX;
	float sideDistY;
	float sideDistZ;

	float deltaDX = abs(1.0f / direction.x);
	float deltaDY = abs(1.0f / direction.y);
	float deltaDZ = abs(1.0f / direction.z);
	float T = -1.0;

	int stepX;
	int stepY;
	int stepZ;

	int hit = 0;
	int side;

	if (direction.x < 0)
	{
		stepX = -1;
		sideDistX = (origin.x - mapX) * deltaDX;
	} 
	
	else 
	{
		stepX = 1;
		sideDistX = (mapX + 1.0 - origin.x) * deltaDX;
	}

	if (direction.y < 0) 
	{
		stepY = -1;
		sideDistY = (origin.y - mapY) * deltaDY;
	} 
	
	else 
	{
		stepY = 1;
		sideDistY = (mapY + 1.0 - origin.y) * deltaDY;
	}

	if (direction.z < 0) 
	{
		stepZ = -1;
		sideDistZ = (origin.z - mapZ) * deltaDZ;
	} 
	
	else 
	{
		stepZ = 1;
		sideDistZ = (mapZ + 1.0 - origin.z) * deltaDZ;
	}

	for (int i = 0; i < 175; i++) 
	{
		if ((mapX >= MapSize.x && stepX > 0) || (mapY >= MapSize.y && stepY > 0) || (mapZ >= MapSize.z && stepZ > 0)) break;
		if ((mapX < 0 && stepX < 0) || (mapY < 0 && stepY < 0) || (mapZ < 0 && stepZ < 0)) break;

		if (sideDistX < sideDistY && sideDistX < sideDistZ) 
		{
			sideDistX += deltaDX;
			mapX += stepX;
			side = 0;
		} 
		
		else if (sideDistY < sideDistX && sideDistY < sideDistZ)
		{
			sideDistY += deltaDY;
			mapY += stepY;
			side = 1;
		} 
		
		else 
		{
			sideDistZ += deltaDZ;
			mapZ += stepZ;
			side = 2;
		}

		float block = GetVoxel(ivec3(mapX, mapY, mapZ));

		if (block != 0) 
		{
			hit = 1;
			blockType = block;

			if (side == 0) 
			{
				T = (mapX - origin.x + (1 - stepX) / 2) / direction.x + t1;
				normal = vec3(1, 0, 0) * -stepX;
			}

			else if (side == 1) 
			{
				T = (mapY - origin.y + (1 - stepY) / 2) / direction.y + t1;
				normal = vec3(0, 1, 0) * -stepY;
			}

			else
			{
				T = (mapZ - origin.z + (1 - stepZ) / 2) / direction.z + t1;
				normal = vec3(0, 0, 1) * -stepZ;
			}

			break;
		}
	}

	return T;
}

*/


// 
bool PointIsInSphere(vec3 point, float radius)
{
	return ((point.x * point.x) + (point.y * point.y) + (point.z * point.z)) < (radius * radius);
}

vec3 RandomPointInUnitSphereRejective() // unused since this is slow asf, it was only ever used as a test, nothing else.
{
	float x, y, z;
	const int accuracy = 10;

	for (int i = 0 ; i < clamp(accuracy, 2, 40); i++)
	{
		x = nextFloat(RNG_SEED, -1.0f, 1.0f);
		y = nextFloat(RNG_SEED, -1.0f, 1.0f);
		z = nextFloat(RNG_SEED, -1.0f, 1.0f);

		if (PointIsInSphere(vec3(x, y, z), 1.0f))
		{
			return vec3(x, y, z);
		}
	}

	return vec3(x, y, z);
}

vec3 CosineSampleHemisphere_2(float u1, float u2)
{
    float r = sqrt(u1);
    float theta = 2.0f * 3.14159265f * u2;
    float x = r * cos(theta);
    float y = r * sin(theta);
    return vec3(x, y, sqrt(max(0.0f, 1 - u1)));
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

/// Direct lighting ///

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

float GetShadowAt(vec3 pos, in vec3 ldir)
{
	vec3 RayDirection = (ldir);
	
	float T = -1.0f;
	 
	vec3 norm;
	float block;

	if (u_APPLY_PLAYER_SHADOW) {
		float ShadowTMIN = -1.0f, ShadowTMAX = -1.0f;
		bool IntersectsPlayer = RayBoxIntersect(u_ViewerPosition + vec3(0.2f, 0.0f, 0.2f), u_ViewerPosition - vec3(0.75f, 1.75f, 0.75f), pos.xyz, RayDirection, ShadowTMIN, ShadowTMAX);
		if (IntersectsPlayer)
		{
			return 1.0f;
		}
	}

	T = VoxelTraversalDF(pos.rgb, RayDirection, norm, block, 150);
	return float(T > 0.0f);
}

vec2 ReprojectShadow(in vec3 world_pos)
{
	vec3 WorldPos = world_pos;

	vec4 ProjectedPosition = u_ShadowProjection * u_ShadowView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
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






// BRDF

// distribution ggx
float GGX_D(float dotNH, float alpha2)
{ 
    float den = (alpha2 - 1.0) * dotNH * dotNH + 1.0;
    return alpha2 / (PI * den * den);
}

vec3 GGX_LIGHT(vec3 n, vec3 v, vec3 l, vec2 RM, vec3 albedo) 
{
    float alpha2 = RM.x * RM.x;
    float dotNL = clamp(dot(n, l), 0.0f, 1.0f);
    float dotNV = clamp(dot(n, v), 0.0f, 1.0f);
    vec3 h = normalize(v + l);
    float dotNH = clamp(dot(n, h), 0.0f, 1.0f);
    float dotLH = clamp(dot(l, h), 0.0f, 1.0f);
    float D = GGX_D(dotNH, alpha2);
    vec3 F0 = RM.y*albedo;
    vec3 F = F0 + (1.0f - F0) * pow(1.9f - dotLH, 5.0f);
    float k = (RM.x + 1.0f) * (RM.x + 1.0f) / 8.0f;
    float G = 1.0f / ((dotNL * (1.0f - k) + k) * (dotNV * (1.0f - k) + k));
    return max(D * F * G, 0.0001f);
}

float InverseSchlick(float f0, float VoH) 
{
    return 1.0 - clamp(f0 + (1.0f - f0) * pow(1.0f - VoH, 5.0f), 0.0f, 1.0f);
}

// hammon diffuse brdf 
float DiffuseHammon(vec3 normal, vec3 viewDir, vec3 lightDir, float roughness)
{
    float nDotL = max(dot(normal, lightDir), 0.0f);

    if (nDotL <= 0.0) 
	{
		return 0.0f; 
	}

	const float OneOverPI = 1.0f/PI;

    float nDotV = max(dot(normal, viewDir), 0.0f);
    float lDotV = max(dot(lightDir, viewDir), 0.0f);
    vec3 halfWay = normalize(viewDir + lightDir);
    float nDotH = max(dot(normal, halfWay), 0.0f);
    float facing = lDotV * 0.5f + 0.5f;
    float singleRough = facing * (0.9f - 0.4f * facing) * ((0.5f + nDotH) * rcp(max(nDotH, 0.02)));
    float singleSmooth = 1.05f * InverseSchlick(0.0f, nDotL) * InverseSchlick(0.0f, max(nDotV, 0.0f));
    float single = clamp(mix(singleSmooth, singleRough, roughness) * rcp(PI), 0.0f, 1.0f);
    float multi = 0.1159f * roughness;
    return clamp((multi + single) * nDotL, 0.0f, 1.0f); // approximate for multi scattering as well
}

vec3 GetBRDF(vec3 normal, vec3 incident, vec3 sunDir, vec3 albedo, vec2 PBR) 
{
    vec3 BRDF;
    bool isMetal = (PBR.y > 0.08f);
    float cosTheta = clamp(dot(normal, sunDir), 0.0f, 1.0f);

    if(!isMetal) 
	{
        return albedo * DiffuseHammon(normal, -incident, sunDir, PBR.x);
    }
	
	else 
	{
       return mix(vec3(1), albedo, float(isMetal)) * GGX_LIGHT(normal, -incident, sunDir, PBR, albedo) * cosTheta;
    }
    
    return BRDF;
}

// End