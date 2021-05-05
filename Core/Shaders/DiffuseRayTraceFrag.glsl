#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

//#define USE_HEMISPHERICAL_DIFFUSE_SCATTERING
#define ANIMATE_NOISE
#define MAX_VOXEL_DIST 10

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler3D u_VoxelData;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_PositionTexture;
uniform samplerCube u_Skymap;
uniform vec2 u_Dimensions;
uniform float u_Time;

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

vec3 g_Normal;
int RNG_SEED = 0;

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

vec3 GetBlockRayColor(in Ray r, out float T, out vec3 out_n, out bool intersection)
{
	float b;

	T = voxel_traversal(r.Origin, r.Direction, out_n, b, 10);

	if (T > 0.0f) 
	{ 
		return vec3(0.7f); 
	} 

	else 
	{	
		//return GetSkyColorAt(r.Direction);
		return vec3(1.0f);
	}
}

vec3 CalculateDiffuse(in vec3 initial_origin, in vec3 initial_normal)
{
	Ray r1, r2;
	float t1, t2;
	vec3 n1, n2;
	bool i1, i2;

	r1.Origin = initial_origin;

	#ifdef USE_HEMISPHERICAL_DIFFUSE_SCATTERING
	r1.Direction = cosWeightedRandomHemisphereDirection(initial_normal);
	#else 
	r1.Direction = initial_normal + RandomPointInUnitSphereRejective();
	#endif

	vec3 col_1 = GetBlockRayColor(r1, t1, n1, i1); 

	r2.Origin = r1.Origin + (r1.Direction * t1);

	#ifdef USE_HEMISPHERICAL_DIFFUSE_SCATTERING
	r2.Direction = n1 + cosWeightedRandomHemisphereDirection(n1);
	#else 
	r2.Direction = n1 + RandomPointInUnitSphereRejective();
	#endif

	vec3 col_2 = GetBlockRayColor(r2, t2, n2, i2); 

	vec3 diff = mix(col_1, col_2, 0.4f);
	diff = pow(diff, vec3(3.5f));

	return vec3(diff);
}

void main()
{
    #ifdef ANIMATE_NOISE
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(u_Time * 1000);
	#else
		RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x);
	#endif

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
		float u = (Pixel.x + nextFloat(RNG_SEED)) / u_Dimensions.x;
		float v = (Pixel.y + nextFloat(RNG_SEED)) / u_Dimensions.y;
		vec2 uv = vec2(u, v);

		vec4 Position = texture(u_PositionTexture, uv); // initial intersection point
		vec3 Normal = texture(u_NormalTexture, uv).rgb;

		o_Color += CalculateDiffuse(Position.xyz, Normal);

	}

	o_Color = o_Color / SPP;
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