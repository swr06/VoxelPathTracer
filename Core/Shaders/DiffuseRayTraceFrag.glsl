#version 330 core


#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

#define USE_COLORED_DIFFUSE // Applies diffuse from the block albedo
#define USE_HEMISPHERICAL_DIFFUSE_SCATTERING 
#define ANIMATE_NOISE // Has to be enabled for temporal filtering to work properly 
#define MAX_VOXEL_DIST 20
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

// Outputs diffuse indirect
// Specular indirect is handled separately and in a higher resolution
layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler3D u_VoxelData;
uniform sampler3D u_DistanceFieldTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_PositionTexture;
uniform samplerCube u_Skymap;

uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlockEmissiveTextures;

uniform sampler2D u_DataTexture;
//uniform sampler2D u_BlueNoiseTexture; // Single 256x256 blue noise texture

uniform vec2 u_Dimensions;
uniform float u_Time;

uniform int u_SPP;

uniform vec3 u_ViewerPosition;

uniform vec4 BLOCK_TEXTURE_DATA[128]; // Albedo tex index, normal tex index, pbr tex index, emissive 

uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;

// Temp
uniform mat4 u_ShadowView;
uniform mat4 u_ShadowProjection;
uniform sampler2D u_ShadowMap;
// Temp

uniform float u_DiffuseLightIntensity = 4.0f;

// Function prototypes
float nextFloat(inout int seed, in float min, in float max);
float nextFloat(inout int seed, in float max);
vec3 cosWeightedRandomHemisphereDirection(const vec3 n); 
float nextFloat(inout int seed); 
int nextInt(inout int seed); 
vec3 GetSkyColorAt(vec3 rd);
vec3 RandomPointInUnitSphereRejective();
vec3 CosineSampleHemisphere(float u1, float u2);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
float GetShadowAt(vec3 pos, in vec3 ldir);
vec2 ReprojectShadow(vec3);
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv, out int NormalIndex);
float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType, in int dist);
bool voxel_traversal(vec3 origin, vec3 direction, inout float block, out vec3 normal, out vec3 world_pos, int dist);


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

// Simplified diffuse brdf
vec3 CalculateDirectionalLight(in vec3 world_pos, in vec3 light_dir, vec3 radiance, in vec3 albedo, in vec3 normal, in float shadow)
{
	vec3 DiffuseBRDF = albedo * max(dot(normal, normalize(light_dir)), 0.0f) * (radiance * 1.5f);
    return DiffuseBRDF * (1.0f - shadow);
} 

vec3 GetDirectLighting(in vec3 world_pos, in int tex_index, in vec2 uv, in vec3 flatnormal)
{
	vec3 SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * (8.0f);
	vec3 NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.8f; 
	vec3 DUSK_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.9f; 

	vec3 LIGHT_COLOR; // The radiance of the light source
	vec3 StrongerLightDirection;
	bool SunStronger = -u_SunDirection.y < 0.01f ? true : false;
	float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 2.9), 0.0f, 1.0f);
    vec3 SunColor = mix(SUN_COLOR, DUSK_COLOR, DuskVisibility);
	LIGHT_COLOR = SunStronger ? SunColor : NIGHT_COLOR;
	StrongerLightDirection = SunStronger ? u_SunDirection : u_MoonDirection;
	vec4 TextureIndexes = BLOCK_TEXTURE_DATA[tex_index].xyzw;
	vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(uv, TextureIndexes.r), 3).rgb;

	float Emmisivity = 0.0f;

	if (TextureIndexes.w >= 0.0f)
	{
		float SampledEmmisivity = texture(u_BlockEmissiveTextures, vec3(uv, TextureIndexes.w)).r;
		Emmisivity = SampledEmmisivity * 20.0f;
	}

	vec3 bias = (flatnormal * 0.045);
	float ShadowAt = GetShadowAt(world_pos + bias, StrongerLightDirection);
	vec3 DirectLighting = CalculateDirectionalLight(world_pos, normalize(StrongerLightDirection), LIGHT_COLOR, Albedo, flatnormal, ShadowAt);
	return (Emmisivity * Albedo) + DirectLighting;
}

vec3 GetBlockRayColor(in Ray r, out vec3 pos, out vec3 out_n)
{
	float b = 0;
	int temp;

	bool Intersect = voxel_traversal(r.Origin, r.Direction, b, out_n, pos, MAX_VOXEL_DIST);
	int tex_ref = clamp(int(floor(b * 255.0f)), 0, 127);

	if (Intersect && b > 0) 
	{ 
		vec2 txc; CalculateUV(pos, out_n, txc, temp);
		return GetDirectLighting(pos, tex_ref, txc, out_n);
	} 

	else 
	{	
		return GetSkyColorAt(r.Direction) * 1.35f;
	}
}

vec4 CalculateDiffuse(in vec3 initial_origin, in vec3 input_normal)
{
	float bias = 0.0f;
	Ray new_ray = Ray(initial_origin + input_normal * bias, cosWeightedRandomHemisphereDirection(input_normal));

	vec3 total_color = vec3(0.0f);;

	vec3 Position;
	vec3 Normal;
	float ao;

	for (int i = 0 ; i < MAX_BOUNCE_LIMIT ; i++)
	{
		

		vec3 tangent_normal;
		total_color += GetBlockRayColor(new_ray, Position, Normal);
		new_ray.Origin = Position + Normal * bias;
		new_ray.Direction = cosWeightedRandomHemisphereDirection(Normal);

		if (i == 0)
		{
			// Calculate ao on first bounce
			float dist = distance(initial_origin, Position);
			ao = 1.0f - float(dist*dist < 0.4f);
		}
	}
	
	total_color = total_color / max(MAX_BOUNCE_LIMIT, 1);
	return vec4(total_color, ao); 
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
	o_Color.w = 0.0f;

	if (InitialTracePosition.w <= 0.0f)
	{
		o_Color.xyz = GetSkyColorAt(normalize(v_RayDirection));
		return;
	}

	vec3 TotalColor = vec3(0.0f);
	vec3 Normal = texture(u_NormalTexture, v_TexCoords).rgb;
	float AccumulatedAO = 0.0f;

	int SPP = clamp(u_SPP, 1, 32);

	for (int s = 0 ; s < SPP ; s++)
	{
		vec4 x = CalculateDiffuse(InitialTracePosition.xyz, Normal);
		TotalColor += x.xyz;
		AccumulatedAO += x.w;
	}

	o_Color.xyz = TotalColor / SPP;
	o_Color.xyz += vec3(0.00525f);
	o_Color.w = AccumulatedAO / SPP;
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

float GetShadowAt(vec3 pos, in vec3 ldir)
{
	vec3 RayDirection = normalize(ldir);
	
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