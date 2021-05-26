#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

#define ALBEDO_TEX_LOD 3 // 512, 256, 128
//#define JITTER_BASED_ON_ROUGHNESS

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_InitialTraceNormalTexture;
uniform sampler2D u_PBRTexture;
uniform sampler2D u_BlueNoiseTexture;

uniform sampler3D u_VoxelData;
uniform samplerCube u_Skymap;

uniform vec4 BLOCK_TEXTURE_DATA[128];
uniform float u_ReflectionTraceRes;

uniform vec2 u_Dimensions;
uniform vec3 u_SunDirection;
uniform vec3 u_StrongerLightDirection;
uniform float u_Time;

uniform int u_GrassBlockProps[10];

uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockPBRTextures;

uniform vec3 u_ViewerPosition;
		
// Function prototypes
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType, in int mdist) ;
float GetVoxel(ivec3 loc);
bool IsInVoxelizationVolume(in vec3 pos);
float GetShadowAt(in vec3 pos, in vec3 ldir);

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

    return clamp(Result, 0.0f, 1000.0f) * clamp((1.0f - Shadow), 0.0f, 1.0f);
}

int BLUE_NOISE_IDX = 0;

float GetBlueNoise()
{
	BLUE_NOISE_IDX++;
	vec2 txc =  vec2(BLUE_NOISE_IDX / 255, mod(BLUE_NOISE_IDX, 255));
	return texelFetch(u_BlueNoiseTexture, ivec2(txc), 0).r;
}

vec3 ImportanceSampleGGX(vec3 N, float roughness)
{
	vec2 Xi;

	Xi.x = GetBlueNoise();
	Xi.y = GetBlueNoise();

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
    if (fragpos.y < 0.24f) { return 0.0f; }

	float elevation = clamp(fragpos.y, 0.0f, 1.0f);
	vec2 uv = fragpos.xz / (1.0f + elevation);

    float star = StableStarField(uv * 700.0f, 0.999);
    
    // Star shimmer
    float rand_val = rand(fragpos.xy);
    star *= (rand_val + sin(u_Time * rand_val) * 1.5f);

	return clamp(star, 0.0f, 100000.0f) * 30.0f;
}

const vec3 ATMOSPHERE_SUN_COLOR = vec3(1.0f * 6.25f, 1.0f * 6.25f, 0.8f * 4.0f);
const vec3 ATMOSPHERE_MOON_COLOR =  vec3(0.7f, 0.7f, 1.25f);

bool GetAtmosphere(inout vec3 atmosphere_color, in vec3 in_ray_dir)
{
    vec3 sun_dir = normalize(u_SunDirection); 
    vec3 moon_dir = vec3(-sun_dir.x, -sun_dir.y, sun_dir.z); 

    vec3 ray_dir = normalize(in_ray_dir);
    vec3 atmosphere = texture(u_Skymap, ray_dir).rgb;
    bool intersect = false;

    if(dot(ray_dir, sun_dir) > 0.9997f)
    {
        atmosphere *= ATMOSPHERE_SUN_COLOR * 3.0f; intersect = true;
    }

    if(dot(ray_dir, moon_dir) > 0.99986f)
    {
        atmosphere *= ATMOSPHERE_MOON_COLOR * 50.0f; intersect = true;
    }

    float star_visibility;
    star_visibility = clamp(exp(-distance(-u_SunDirection.y, 1.8555f)), 0.0f, 1.0f);
    vec3 stars = vec3(stars(vec3(in_ray_dir)) * star_visibility);
    stars = clamp(stars, 0.0f, 1.3f);

    atmosphere += stars;

    atmosphere_color = atmosphere;

    return intersect;
}

const vec3 SUN_COLOR = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * 7.4f;
const vec3 SUN_AMBIENT = (vec3(120.0f, 172.0f, 255.0f) / 255.0f) * 1.01f;

const vec3 NIGHT_COLOR  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 2.1f; 
const vec3 NIGHT_AMBIENT  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.76f; 

const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

void main()
{
	// Start ray at sampled position, use normalized normal (already in tangent space) as direction, trace and get the albedo color at.
	
	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * 800 * int(floor(u_Time * 100));
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
	BLUE_NOISE_IDX = int(floor(RNG_SEED));
	BLUE_NOISE_IDX = BLUE_NOISE_IDX % (255 * 255);

	vec2 Pixel;
	Pixel.x = v_TexCoords.x * u_Dimensions.x;
	Pixel.y = v_TexCoords.y * u_Dimensions.y;

	const int SPP = 2;
	int total_hits = 0;
	vec3 TotalColor = vec3(0.0f);

	for (int s = 0 ; s < SPP ; s++)
	{
		vec2 suv;

		float u = (Pixel.x + GetBlueNoise()) / u_Dimensions.x;
		float v = (Pixel.y + GetBlueNoise()) / u_Dimensions.y;

		suv = vec2(u, v);

		vec4 SampledWorldPosition = texture(u_PositionTexture, suv); // initial intersection point
		vec3 InitialTraceNormal = texture(u_InitialTraceNormalTexture, suv).rgb;
		vec4 data = texture(u_PBRTexture, suv);

		vec2 iUV;
		CalculateUV(SampledWorldPosition.xyz, InitialTraceNormal, iUV);
		
		vec4 PBRMap = texture(u_BlockPBRTextures, vec3(iUV, data.z)).rgba;
		float RoughnessAt = PBRMap.r;
		float MetalnessAt = PBRMap.g;

		vec3 ReflectionNormal = ImportanceSampleGGX(InitialTraceNormal, RoughnessAt * 0.75f);

		vec3 I = normalize(SampledWorldPosition.xyz - u_ViewerPosition);
		vec3 R = normalize(reflect(I, ReflectionNormal));
		vec3 Normal;
		float Blocktype;

		float T = voxel_traversal(SampledWorldPosition.xyz, R, Normal, Blocktype, 100);
		vec3 HitPosition = SampledWorldPosition.xyz + (R * T);

		vec2 UV; 
		vec3 Tangent, Bitangent;
		CalculateVectors(HitPosition, Normal, Tangent, Bitangent, UV); UV.y = 1.0f - UV.y;

		if (T > 0.0f)
		{
			int reference_id = clamp(int(floor(Blocktype * 255.0f)), 0, 127);
			vec4 texture_ids = BLOCK_TEXTURE_DATA[reference_id];

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
			TBN = mat3(normalize(Tangent), normalize(Bitangent), normalize(Normal));

			vec3 Albedo = textureLod(u_BlockAlbedoTextures, vec3(UV,texture_ids.x), ALBEDO_TEX_LOD).rgb;
			bool SunStronger = u_StrongerLightDirection == u_SunDirection;
			vec3 Radiance = SunStronger ? SUN_COLOR : NIGHT_COLOR; Radiance *= 3.56f;
			vec3 Ambient = SunStronger ? SUN_AMBIENT : NIGHT_AMBIENT;
			Ambient = (Albedo * Ambient);
				
			vec4 SampledPBR = textureLod(u_BlockPBRTextures, vec3(UV, texture_ids.z), 2).rgba;
			float AO = pow(SampledPBR.w, 2.0f);

			vec3 NormalMapped = TBN * (textureLod(u_BlockNormalTextures, vec3(UV,texture_ids.y), 2).rgb * 2.0f - 1.0f);
			vec3 DirectLighting = (Ambient * 0.6f) + 
									CalculateDirectionalLight(HitPosition, 
																u_StrongerLightDirection, 
																Radiance, 
																Albedo, 
																NormalMapped, 
																SampledPBR.xyz,
																GetShadowAt(HitPosition, u_StrongerLightDirection) * 0.9);
			
			vec3 Computed;
			Computed = DirectLighting;
			Computed *= AO;
			TotalColor += Computed;

		}

		else
		{
			bool BodyIntersect;
			vec3 AtmosphereColor;

			vec3 AtmosphereRayDir = R;
			AtmosphereRayDir.y = clamp(R.y, 0.1f, 1.5f);

			BodyIntersect = GetAtmosphere(AtmosphereColor, AtmosphereRayDir);
			TotalColor += AtmosphereColor * 1.6f;
		}


		total_hits++;
	}

	o_Color = (TotalColor / float(total_hits));
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

			int reference_id = clamp(int(floor(block * 255.0f)), 0, 127);
			bool transparent = BLOCK_TEXTURE_DATA[reference_id].a > 0.5f;

			if (transparent)
			{
				vec3 hit_position = orig + (direction * T);
				vec2 uv;

				CalculateUV(hit_position, normal, uv); uv.y = 1.0f - uv.y;

				if (textureLod(u_BlockAlbedoTextures, vec3(uv, BLOCK_TEXTURE_DATA[reference_id].x), 1).a < 0.05f)
				{
					T = -1.0f;
					continue;
				}
			}

			break;
		}
	}

	return T;
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

void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );

	if (normal == Normals[0])
    {
        uv = vec2(fract(world_pos.xy));
    }

    else if (normal == Normals[1])
    {
        uv = vec2(fract(world_pos.xy));
    }

    else if (normal == Normals[2])
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (normal == Normals[3])
    {
        uv = vec2(fract(world_pos.xz));
    }
	
    else if (normal == Normals[4])
    {
        uv = vec2(fract(world_pos.zy));
    }
    

    else if (normal == Normals[5])
    {
        uv = vec2(fract(world_pos.zy));
    }
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

float GetShadowAt(in vec3 pos, in vec3 ldir)
{
	vec3 RayDirection = normalize(ldir);
	
	float T = -1.0f;
	 
	vec3 norm;
	float block;

	float ShadowTMIN = -1.0f, ShadowTMAX = -1.0f;
	bool PlayerIntersect = RayBoxIntersect(u_ViewerPosition + vec3(0.2f, 0.0f, 0.2f), u_ViewerPosition - vec3(0.75f, 1.75f, 0.75f), pos.xyz, RayDirection, ShadowTMIN, ShadowTMAX);
	if (PlayerIntersect) { return 1.0f; }

	T = voxel_traversal(pos.rgb, RayDirection, norm, block, 40);

	if (T > 0.0f) 
	{ 
		return 1.0f; 
	}
	
	else 
	{ 
		return 0.0f;
	}
}