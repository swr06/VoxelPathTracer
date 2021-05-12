#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

#define USE_COLORED_DIFFUSE // Applies diffuse from the block albedo
#define USE_HEMISPHERICAL_DIFFUSE_SCATTERING 
//#define USE_BAYER_PIXEL_DITHER
#define ANIMATE_NOISE // Has to be enabled for temporal filtering to work properly 
#define MAX_VOXEL_DIST 20 // Increase this when DAGs are implemented 
#define NORMAL_MAP_LOD 3 // 512, 256, 128, 64, 32, 16, 8, 4, 2
#define ALBEDO_MAP_LOD 4 // 512, 256, 128, 64, 32, 16, 8, 4, 2
#define MAX_BOUNCE_LIMIT 4

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

uniform sampler2D u_DataTexture;

uniform vec2 u_Dimensions;
uniform float u_Time;

uniform vec3 BLOCK_TEXTURE_DATA[128];


// Function prototypes
float nextFloat(inout int seed, in float min, in float max);
float nextFloat(inout int seed, in float max);
vec3 cosWeightedRandomHemisphereDirection(const vec3 n); 
float nextFloat(inout int seed); 
int nextInt(inout int seed); 
vec3 GetSkyColorAt(vec3 rd);
float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType, in int);
float ProjectToCube(vec3 ro, vec3 rd) ;
bool VoxelExists(in vec3 loc);
float GetVoxel(ivec3 loc);
bool IsInVoxelizationVolume(in vec3 pos);
vec3 RandomPointInUnitSphereRejective();
vec3 CosineSampleHemisphere(float u1, float u2);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);

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

vec3 GetBlockRayColor(in Ray r, out float T, out vec3 out_n, out bool intersection, 
					  out int tex_ref, out vec3 tangent, out vec3 bitangent, out vec2 txc)
{
	float b;

	T = voxel_traversal(r.Origin, r.Direction, out_n, b, MAX_VOXEL_DIST);
	tex_ref = clamp(int(floor(b * 255.0f)), 0, 127);

	if (T > 0.0f) 
	{ 
		CalculateVectors(r.Origin + (r.Direction * T), out_n, tangent, bitangent, txc);

		#ifdef USE_COLORED_DIFFUSE
		int albedo_tex_index = int(floor(BLOCK_TEXTURE_DATA[tex_ref].r));
		return textureLod(u_BlockAlbedoTextures, vec3(txc, albedo_tex_index), ALBEDO_MAP_LOD).rgb * 0.96f;
		#else 
		return vec3(0.2f, 0.2f, 0.23f);
		#endif
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
		total_color = GetBlockRayColor(new_ray, T, N, intersection, tex_ref, Tangent, Bitangent, TexCoord);
		tbn =  mat3(normalize(Tangent), normalize(Bitangent), normalize(N));

		vec3 tangent_normal = tbn * (textureLod(u_BlockNormalTextures, vec3(TexCoord, BLOCK_TEXTURE_DATA[tex_ref].g), NORMAL_MAP_LOD).rgb * 2.0f - 1.0f);
		new_ray.Origin = new_ray.Origin + (new_ray.Direction * T);
		new_ray.Direction = cosWeightedRandomHemisphereDirection(tangent_normal);
	}
	
	total_color = pow((total_color), vec3(MAX_BOUNCE_LIMIT - 1));
	//total_color = total_color / MAX_BOUNCE_LIMIT;

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
	#ifdef USE_BAYER_PIXEL_DITHER
		float u = (Pixel.x + Bayer16(vec2(Pixel.x + s, Pixel.y * s))) / u_Dimensions.x;
		float v = (Pixel.y + Bayer16(vec2(Pixel.x + s, Pixel.y * s))) / u_Dimensions.y;
	#else 
		float u = (Pixel.x + nextFloat(RNG_SEED)) / u_Dimensions.x;
		float v = (Pixel.y + nextFloat(RNG_SEED)) / u_Dimensions.y;
	#endif

		vec2 uv = vec2(u, v);

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
  	vec2 r = vec2(nextFloat(RNG_SEED), nextFloat(RNG_SEED));
    
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
