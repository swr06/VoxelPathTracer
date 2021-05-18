#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

#define USE_COLORED_DIFFUSE // Applies diffuse from the block albedo
#define USE_HEMISPHERICAL_DIFFUSE_SCATTERING 
//#define USE_BAYER_PIXEL_DITHER
#define ANIMATE_NOISE // Has to be enabled for temporal filtering to work properly 
#define MAX_VOXEL_DIST 16 
#define MAX_SHADOW_TRACE_DIST 150 // Needs to be high
#define NORMAL_MAP_LOD 3 // 512, 256, 128, 64, 32, 16, 8, 4, 2
#define ALBEDO_MAP_LOD 4 // 512, 256, 128, 64, 32, 16, 8, 4, 2
#define MAX_BOUNCE_LIMIT 2
#define APPLY_PLAYER_SHADOW

// Bayer matrix, used for testing dithering
#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler3D u_VoxelData;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_PositionTexture;
uniform samplerCube u_Skymap;

uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlockEmissiveTextures;

uniform sampler2D u_DataTexture;
uniform sampler2D u_BlueNoiseTexture; // Single 256x256 blue noise texture

uniform vec2 u_Dimensions;
uniform float u_Time;

uniform vec3 u_ViewerPosition;

uniform vec4 BLOCK_TEXTURE_DATA[128]; // Albedo tex index, normal tex index, pbr tex index, emissive 

uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;

// Temp
uniform mat4 u_ShadowView;
uniform mat4 u_ShadowProjection;
uniform sampler2D u_ShadowMap;
// Temp

// Function prototypes
float nextFloat(inout int seed, in float min, in float max);
float nextFloat(inout int seed, in float max);
vec3 cosWeightedRandomHemisphereDirection(const vec3 n); 
float nextFloat(inout int seed); 
int nextInt(inout int seed); 
vec3 GetSkyColorAt(vec3 rd);
float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType, in int mdist);
float ProjectToCube(vec3 ro, vec3 rd) ;
bool VoxelExists(in vec3 loc);
float GetVoxel(ivec3 loc);
bool IsInVoxelizationVolume(in vec3 pos);
vec3 RandomPointInUnitSphereRejective();
vec3 CosineSampleHemisphere(float u1, float u2);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
float GetShadowAt(in vec3 pos, in vec3 ldir);
vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 albedo, vec3 normal, vec3 pbr, float shadow);
vec2 ReprojectShadow(vec3);
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv, out int NormalIndex);

// Globals
vec3 g_Normal;
int RNG_SEED = 0;
int BLUE_NOISE_IDX = 0;

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

float GetBlueNoise()
{
	BLUE_NOISE_IDX++;
	vec2 txc =  vec2(BLUE_NOISE_IDX / 256, mod(BLUE_NOISE_IDX, 256));
	return texelFetch(u_BlueNoiseTexture, ivec2(txc), 0).r;
}

vec3 GetDirectLighting(in vec3 world_pos, in int tex_index, in vec3 normal, in vec2 uv)
{
	vec3 LIGHT_COLOR; // The radiance of the light source
	vec3 StrongerLightDirection;
	bool SunStronger = -u_SunDirection.y < 0.01f ? true : false;
	LIGHT_COLOR = SunStronger ? vec3(1.0f) * 18.0f : vec3(1.0f) * 6.0f;
	StrongerLightDirection = SunStronger ? u_SunDirection : u_MoonDirection;

	vec4 TextureIndexes = BLOCK_TEXTURE_DATA[tex_index].xyzw;
	vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(uv, TextureIndexes.r), 2).rgb;
	vec3 PBR = textureLod(u_BlockPBRTextures, vec3(uv, TextureIndexes.b), 2).rgb;

	float Emmisivity = 0.0f;

	if (TextureIndexes.w >= 0.0f)
	{
		float SampledEmmisivity = texture(u_BlockEmissiveTextures, vec3(uv, TextureIndexes.w)).r;
		Emmisivity = SampledEmmisivity * 65.0f;
	}

	float ShadowAt = GetShadowAt(world_pos, StrongerLightDirection);
	vec3 DirectLighting = CalculateDirectionalLight(world_pos, normalize(StrongerLightDirection), LIGHT_COLOR, Albedo, normal, PBR, ShadowAt);
	return (Emmisivity * Albedo) + DirectLighting;
}

vec3 GetBlockRayColor(in Ray r, out float T, out vec3 out_n, out bool intersection, 
					  out int tex_ref, out vec3 tangent, out vec3 bitangent, out vec2 txc, out vec3 tangent_normal)
{
	float b;

	T = voxel_traversal(r.Origin, r.Direction, out_n, b, MAX_VOXEL_DIST);
	tex_ref = clamp(int(floor(b * 255.0f)), 0, 127);

	if (T > 0.0f) 
	{ 
		CalculateVectors(r.Origin + (r.Direction * T), out_n, tangent, bitangent, txc);

		mat3 tbn =  mat3(normalize(tangent), normalize(bitangent), normalize(out_n));
		tangent_normal = tbn * (textureLod(u_BlockNormalTextures, vec3(txc, BLOCK_TEXTURE_DATA[tex_ref].g), NORMAL_MAP_LOD).rgb * 2.0f - 1.0f);
		
		return GetDirectLighting(r.Origin + (r.Direction * T), tex_ref, tangent_normal, txc);
	} 

	else 
	{	
		return GetSkyColorAt(r.Direction);
	}
}

vec3 CalculateDiffuse(in vec3 initial_origin, in vec3 input_normal)
{
	vec3 initial_normal;
	vec2 initial_uv;
	vec3 initial_tan;
	vec3 initial_bitan;
	mat3 initial_tbn;
	int initial_idx;

	CalculateVectors(initial_origin, input_normal, initial_tan, initial_bitan, initial_uv);
	initial_tbn = mat3(normalize(initial_tan), normalize(initial_bitan), normalize(input_normal));
	initial_idx = int(floor(texture(u_DataTexture, v_TexCoords).g));
	initial_normal = initial_tbn * (textureLod(u_BlockNormalTextures, vec3(initial_uv, initial_idx), NORMAL_MAP_LOD).rgb * 2.0f - 1.0f);

	Ray new_ray = Ray(initial_origin, cosWeightedRandomHemisphereDirection(initial_normal));

	vec3 total_color = vec3(0.0f);;

	float T;
	vec3 N;
	bool intersection;
	int tex_ref;
	vec3 Tangent, Bitangent;
	vec2 TexCoord;
	mat3 tbn;

	for (int i = 0 ; i < MAX_BOUNCE_LIMIT ; i++)
	{
		vec3 tangent_normal;
		total_color += GetBlockRayColor(new_ray, T, N, intersection, tex_ref, Tangent, Bitangent, TexCoord, tangent_normal);
		new_ray.Origin = new_ray.Origin + (new_ray.Direction * T);
		new_ray.Direction = cosWeightedRandomHemisphereDirection(tangent_normal);
	}
	
	//total_color = pow((total_color), vec3(3.0f)) / MAX_BOUNCE_LIMIT;
	total_color = total_color / max(MAX_BOUNCE_LIMIT - 0.25f, 0.66f);

	return total_color;
}

void main()
{
    #ifdef ANIMATE_NOISE
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(u_Time * 1000);
	#else
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x);
	#endif

	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;

	BLUE_NOISE_IDX += int(floor(RNG_SEED));
	BLUE_NOISE_IDX = BLUE_NOISE_IDX % (255 * 255);

	vec4 InitialTracePosition = texture(u_PositionTexture, v_TexCoords).rgba;

	if (InitialTracePosition.w == -1.0f)
	{
		o_Color = GetSkyColorAt(normalize(v_RayDirection));
		return;
	}

    vec3 Normals[6] = vec3[](
	                  vec3(0.0f, 1.0f, 0.0f),
                      vec3(0.0f, -1.0f, 0.0f),
                      vec3(0.0f, 0.0f, 1.0f),
                      vec3(0.0f, 0.0f, -1.0f),
                      vec3(-1.0f, 0.0f, 0.0f),
                      vec3(1.0f, 0.0f, 0.0f)
    );

	vec2 Pixel;
	Pixel.x = v_TexCoords.x * u_Dimensions.x;
	Pixel.y = v_TexCoords.y * u_Dimensions.y;

	const int SPP = 4;

	for (int s = 0 ; s < SPP ; s++)
	{
		vec2 uv;

		if (s == 0)
		{
			float u = (Pixel.x) / u_Dimensions.x;
			float v = (Pixel.y) / u_Dimensions.y;

			uv = vec2(u, v);
		}

		else 
		{
			float u = (Pixel.x + GetBlueNoise()) / u_Dimensions.x;
			float v = (Pixel.y + GetBlueNoise()) / u_Dimensions.y;

			uv = vec2(u, v);
		}


		vec4 Position = texture(u_PositionTexture, uv); // initial intersection point
		vec3 Normal = texture(u_NormalTexture, uv).rgb;

		o_Color += CalculateDiffuse(Position.xyz, Normal);

	}

	o_Color = o_Color / SPP;
	o_Color = clamp(o_Color, 0.02f, 100000.0f);
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
  	vec2 r = vec2(GetBlueNoise(), GetBlueNoise());
    
	vec3  uu = normalize(cross(n, vec3(0.0,1.0,1.0)));
	vec3  vv = cross(uu, n);
	float ra = sqrt(r.y);
	float rx = ra * cos(6.2831 * r.x); 
	float ry = ra * sin(6.2831 * r.x);
	float rz = sqrt(1.0 - r.y);
	vec3  rr = vec3(rx * uu + ry * vv + rz * n );
    
    return normalize(rr);
}

vec3 GetSkyColorAt(vec3 rd) 
{
    return texture(u_Skymap, (rd)).rgb;
}

bool IsInVoxelizationVolume(in vec3 pos)
{
    if (pos.x < 0.0f || pos.y < 0.0f || pos.z < 0.0f || 
        pos.x > float(WORLD_SIZE_X) || pos.y > float(WORLD_SIZE_Y) || pos.z > float(WORLD_SIZE_Z))
    {
        return false;    
    }   

    return true;
}

float GetVoxel(ivec3 loc)
{
    if (IsInVoxelizationVolume(loc))
    {
         return texelFetch(u_VoxelData, loc, 0).r;
    }
    
    return 0.0f;
}

bool VoxelExists(in vec3 loc)
{
    if (GetVoxel(ivec3(loc)) > 0.0f) 
    {
        return true;
    }

    return false;
}

float ProjectToCube(vec3 ro, vec3 rd) 
{	
	const vec3 MapSize = vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);
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

float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType, in int mdist) 
{
	const vec3 MapSize = vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);

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

	for (int i = 0; i < mdist; i++) 
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

			//int reference_id = clamp(int(floor(block * 255.0f)), 0, 127);
			//bool transparent = BLOCK_TEXTURE_DATA[reference_id].a > 0.5f;
			//
			//if (transparent)
			//{
			//	vec3 hit_position = orig + (direction * T);
			//	int temp_idx; 
			//	vec2 uv;
			//
			//	CalculateUV(hit_position, normal, uv, temp_idx); uv.y = 1.0f - uv.y;
			//
			//	if (textureLod(u_BlockAlbedoTextures, vec3(uv, BLOCK_TEXTURE_DATA[reference_id].x), 3).a < 0.1f)
			//	{
			//		T = -1.0f;
			//		continue;
			//	}
			//}

			break;
		}
	}

	return T;
}

bool PointIsInSphere(vec3 point, float radius)
{
	return ((point.x * point.x) + (point.y * point.y) + (point.z * point.z)) < (radius * radius);
}

vec3 RandomPointInUnitSphereRejective()
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

vec3 CosineSampleHemisphere(float u1, float u2)
{
    float r = sqrt(u1);
    float theta = 2 * 3.14159265 * u2;
 
    float x = r * cos(theta);
    float y = r * sin(theta);
 
    return vec3(x, y, sqrt(max(0.0f, 1 - u1)));
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

vec2 ReprojectShadow(in vec3 world_pos)
{
	vec3 WorldPos = world_pos;

	vec4 ProjectedPosition = u_ShadowProjection * u_ShadowView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

float GetShadowAt(in vec3 pos, in vec3 ldir)
{
	//vec2 ReprojectedShadow = ReprojectShadow(pos);
	//
	//// Check if in screen space bounds
	//if (ReprojectedShadow.x > 0.0f && ReprojectedShadow.x < 1.0f && ReprojectedShadow.y > 0.0f && ReprojectedShadow.y < 1.0f)
	//{
	//	return texture(u_ShadowMap, ReprojectedShadow).r;
	//}
	
	vec3 RayDirection = normalize(ldir - (ldir * 0.1f));
	
	float T = -1.0f;
	 
	vec3 norm;
	float block;

	#ifdef APPLY_PLAYER_SHADOW
		float ShadowTMIN = -1.0f, ShadowTMAX = -1.0f;
		bool PlayerIntersect = RayBoxIntersect(u_ViewerPosition + vec3(0.2f, 0.0f, 0.2f), u_ViewerPosition - vec3(0.75f, 1.75f, 0.75f), pos.xyz, RayDirection, ShadowTMIN, ShadowTMAX);
		if (PlayerIntersect) { return 1.0f; }
	#endif

	T = voxel_traversal(pos.rgb, RayDirection, norm, block, MAX_SHADOW_TRACE_DIST);

	if (T > 0.0f) 
	{ 
		return 1.0f; 
	}
	
	else 
	{ 
		return 0.0f;
	}
}

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH * NdotH;

    float nom   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return nom / max(denom, 0.001); 
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r * r) / 8.0;

    float nom   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 albedo, vec3 normal, vec3 pbr, float shadow)
{
    float Shadow = min(shadow, 1.0f);

	vec3 V = normalize(u_ViewerPosition - world_pos);
    vec3 L = normalize(light_dir);
    vec3 H = normalize(V + L);

    float Roughness = pbr.r;
    float Metalness = pbr.g;

    float NDF = DistributionGGX(normal, H, Roughness);   
    float G = GeometrySmith(normal, V, L, Roughness);      
    vec3 F = fresnelSchlick(clamp(dot(H, V), 0.0, 1.0), vec3(0.04));
       
    vec3 nominator = NDF * G * F; 
    float denominator = 4.0f * max(dot(normal, V), 0.0) * max(dot(normal, L), 0.0);
    vec3 specular = nominator / max(denominator, 0.001f);
    
    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - Metalness;	
    kD = clamp(kD, 0.0f, 1.0f);
    specular = clamp(specular, 0.0f, 1.0f);

    float NdotL = max(dot(normal, L), 0.0);
	vec3 Result = (kD * albedo / PI + (specular)) * radiance * NdotL;

    return clamp(Result, 0.0f, 2.5f) * clamp((1.0f - Shadow), 0.0f, 1.0f);
}


void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv, out int NormalIndex)
{
	const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
	const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
	const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
	const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
	const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
	const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

    if (normal == NORMAL_TOP)
    {
        uv = vec2(fract(world_pos.xz));
		NormalIndex = 0;
    }

    else if (normal == NORMAL_BOTTOM)
    {
        uv = vec2(fract(world_pos.xz));
		NormalIndex = 1;
    }

    else if (normal == NORMAL_RIGHT)
    {
        uv = vec2(fract(world_pos.zy));
		NormalIndex = 2;
    }

    else if (normal == NORMAL_LEFT)
    {
        uv = vec2(fract(world_pos.zy));
		NormalIndex = 3;
    }
    
    else if (normal == NORMAL_FRONT)
    {
        uv = vec2(fract(world_pos.xy));
		NormalIndex = 4;
    }

     else if (normal == NORMAL_BACK)
    {
        uv = vec2(fract(world_pos.xy));
		NormalIndex = 5;
    }
}