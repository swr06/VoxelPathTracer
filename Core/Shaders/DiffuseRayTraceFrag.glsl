#version 450 core

// indirect diffuse raytracing


#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

#define USE_COLORED_DIFFUSE // Applies diffuse from the block albedo
#define USE_HEMISPHERICAL_DIFFUSE_SCATTERING 
#define ANIMATE_NOISE // Has to be enabled for temporal filtering to work properly 
#define MAX_VOXEL_DIST 30
#define MAX_BOUNCE_LIMIT 2
//#define APPLY_PLAYER_SHADOW

// Bayer matrix, used for testing dithering, unused 
#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

// Outputs diffuse indirect
// Specular indirect is handled separately and in a higher resolution

layout (location = 0) out vec4 o_SH; // Projected radiance onto the first 2 spherical harmonics.
layout (location = 1) out vec2 o_CoCg; // Stores the radiance color data in YCoCg
layout (location = 2) out float o_Utility; // Using the first SH band to get the luminance is slightly erraneous so I store this.
layout (location = 3) out float o_AO; // Basic VXAO

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
uniform float u_Time;

uniform int u_SPP;
uniform int u_CheckerSPP;
uniform int u_CurrentFrame;
uniform int u_CurrentFrameMod512;
uniform int u_CurrentFrameMod128;
uniform bool u_UseBlueNoise;

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
	int sobol_256spp_256d[256*256];
	int scramblingTile[128*128*8];
	int rankingTile[128*128*8];
};

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

// Simplified diffuse brdf
vec3 CalculateDirectionalLight(in vec3 world_pos, vec3 radiance, in vec3 albedo, in float shadow, float NDotL, vec3 PBR, vec3 N, vec3 rD, vec3 sdir)
{
	if (!CAUSTICS) {
		return albedo * DiffuseHammon(N, -rD, sdir, PBR.x) * (radiance*3.5f) * (1.0f - shadow) * PI;
	}

	else {
		return (GetBRDF(N, -rD, sdir, albedo, PBR.xy) * (radiance*4.0f)) * (1.0f - shadow);
	}
} 

vec3 LIGHT_COLOR;
vec3 StrongerLightDirection;

vec3 GetDirectLighting(in vec3 world_pos, in int tex_index, in vec2 uv, in vec3 flatnormal, vec3 rD)
{
	vec2 TextureIndexes = vec2(
		float(BlockAlbedoData[tex_index]),
		float(BlockEmissiveData[tex_index])
	);

	vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(uv, TextureIndexes.r), 3.0f).rgb; // 512, 256, 128, 64, 32, 16
	vec3 PBR = textureLod(u_BlockPBRTextures, vec3(uv, TextureIndexes.r), 2.0f).rgb; // 512, 256, 128, 64, 32, 16

	float Emmisivity = 0.0f;

	if (TextureIndexes.y >= 0.0f)
	{
		float SampledEmmisivity = texture(u_BlockEmissiveTextures, vec3(uv, TextureIndexes.y)).r;
		Emmisivity = SampledEmmisivity * 20.0f * u_DiffuseLightIntensity;
	}

	float NDotL = max(dot(flatnormal, StrongerLightDirection), 0.0f);
	vec3 bias = (flatnormal * 0.045);
	float ShadowAt = NDotL < 0.001f ? 0.0f : GetShadowAt(world_pos + bias, StrongerLightDirection);
	vec3 DirectLighting = CalculateDirectionalLight(world_pos, LIGHT_COLOR, Albedo, ShadowAt, NDotL, PBR, flatnormal, rD, StrongerLightDirection);
	return (Emmisivity * Albedo) + DirectLighting;
}

vec3 GetBlockRayColor(in Ray r, out vec3 pos, out vec3 out_n)
{
	float b = 0;

	//bool Intersect = voxel_traversal(r.Origin, r.Direction, b, out_n, pos, MAX_VOXEL_DIST);

	// float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType, in int dist);
	float T = VoxelTraversalDF(r.Origin, r.Direction, out_n, b, MAX_VOXEL_DIST);
	int tex_ref = clamp(int(floor(b * 255.0f)), 0, 127);
	bool Intersect = T > 0.0f;
	pos = r.Origin + (r.Direction * T);

	if (Intersect && b > 0) 
	{ 
		vec2 txc; CalculateUV(pos, out_n, txc);
		return GetDirectLighting(pos, tex_ref, txc, out_n, r.Direction);
	} 

	else 
	{	
		float x = mix(3.0f, 1.75f, u_SunVisibility);
		x = clamp(x,0.0f,5.0f);
		return GetSkyColorAt(r.Direction) * x;
	}
}

vec4 CalculateDiffuse(in vec3 initial_origin, in vec3 input_normal, out vec3 dir)
{
	float bias = 0.045f;
	Ray new_ray = Ray(initial_origin + input_normal * bias, cosWeightedRandomHemisphereDirection(input_normal));

	vec3 total_color = vec3(0.0f);;

	vec3 Position;
	vec3 Normal;
	float ao = 1.0f;

	for (int i = 0 ; i < MAX_BOUNCE_LIMIT ; i++)
	{
		vec3 tangent_normal;
		total_color += GetBlockRayColor(new_ray, Position, Normal);
		float T = distance(initial_origin, Position);

		if (i == 0)
		{
			if (T < 3.5f && T > 0.0f) 
			{
				// Calculate ao on first bounce
				ao = 1.0f - float(T*T<0.225f);
				//ao = T / 3.5f;
			}

			// store sh direction
			dir = new_ray.Direction;
		}

		new_ray.Origin = Position + Normal * bias;
		new_ray.Direction = cosWeightedRandomHemisphereDirection(Normal);
	}
	
	total_color = total_color / max(MAX_BOUNCE_LIMIT, 1);
	return vec4(total_color, ao); 
}

// basic fract(sin) pseudo random number generator
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
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

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
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
	vec3 SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * (10.0f);
	vec3 NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 1.5f; 
	vec3 DUSK_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.9f; 

	bool SunStronger = -u_SunDirection.y < 0.01f ? true : false;
	float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 2.9), 0.0f, 1.0f);
    vec3 SunColor = mix(SUN_COLOR, DUSK_COLOR, DuskVisibility);
	LIGHT_COLOR = SunStronger ? SunColor : NIGHT_COLOR;
	StrongerLightDirection = SunStronger ? u_SunDirection : u_MoonDirection;

	
    #ifdef ANIMATE_NOISE
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(fract(u_Time) * 1000);
	#else
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x);
	#endif

	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(u_Time) * 100.0f;

	vec4 Position = GetPositionAt(u_PositionTexture, v_TexCoords);
	vec3 Normal = GetNormalFromID(texture(u_NormalTexture, v_TexCoords).x);

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

	vec4 sh_data1 = vec4(0.0f);
	vec2 color_data = vec2(0.0f);
	vec3 radiance = vec3(0.0f);

	for (int s = 0 ; s < SPP ; s++)
	{
		vec3 d = vec3(0.0f);
		vec4 x = CalculateDiffuse(Position.xyz, Normal, d);
		radiance += x.xyz;
		AccumulatedAO += x.w;
		
		float SH[6] = IrridianceToSH(x.xyz, d);
		sh_data1 += vec4(SH[0], SH[1], SH[2], SH[3]);
		color_data += vec2(SH[4], SH[5]);
	}

	AccumulatedAO /= SPP;
	sh_data1 /= SPP;
	color_data /= SPP;
	radiance /= SPP;
	o_SH = sh_data1;
	o_CoCg.xy = color_data;
	o_Utility = GetLuminance(radiance);
	o_AO = AccumulatedAO;
}

vec3 lerp(vec3 v1, vec3 v2, float t)
{
	return (1.0f - t) * v1 + t * v2;
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
// thanks to telo for the help

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

float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType) 
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

vec3 RandomPointInUnitSphereRejective() // unused since this is slow asf
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

	#ifdef APPLY_PLAYER_SHADOW
		float ShadowTMIN = -1.0f, ShadowTMAX = -1.0f;
		bool PlayerIntersect = RayBoxIntersect(u_ViewerPosition + vec3(0.2f, 0.0f, 0.2f), u_ViewerPosition - vec3(0.75f, 1.75f, 0.75f), pos.xyz, RayDirection, ShadowTMIN, ShadowTMAX);
		if (PlayerIntersect) { return 1.0f; }
	#endif

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
float GGX_D(float dotNH, float alpha2) { // distribution
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