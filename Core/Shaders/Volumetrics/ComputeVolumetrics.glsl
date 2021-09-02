#version 430 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

// Bayer dither
#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))


layout (location = 0) out vec3 o_Volumetrics;

in vec2 v_TexCoords;

layout(r16ui, binding = 0) uniform uimage3D u_ParticipatingMedia;
uniform sampler2D u_BlueNoise;
uniform sampler2D u_LinearDepthTexture;

uniform vec3 u_ViewerPosition;
uniform vec2 u_Dimensions;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform int u_LightCount;

layout (std430, binding = 2) buffer SSBO_BlockAverageData
{
    vec4 BlockAverageColorData[128];
};

layout (std430, binding = 4) buffer SSBO_LightData
{
    vec4 LightLocations[1024];
};

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
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

uint GetVoxel(ivec3 loc)
{
    if (IsInVolume(loc))
    {
        return imageLoad(u_ParticipatingMedia, loc).x;
    }
    
    return 0;
}

float saturate(float x) {
	return clamp(x, 1e-5f, 1.0f);
}

// Used to create variation
float mod289(float x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 mod289(vec4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 perm(vec4 x){return mod289(((x * 34.0) + 1.0) * x);}

float noise(vec3 p)
{
    vec3 a = floor(p);
    vec3 d = p - a;
    d = d * d * (3.0 - 2.0 * d);
    vec4 b = a.xxyy + vec4(0.0, 1.0, 0.0, 1.0);
    vec4 k1 = perm(b.xyxy);
    vec4 k2 = perm(k1.xyxy + b.zzww);
    vec4 c = k2 + a.zzzz;
    vec4 k3 = perm(c);
    vec4 k4 = perm(c + 1.0);
    vec4 o1 = fract(k3 * (1.0 / 41.0));
    vec4 o2 = fract(k4 * (1.0 / 41.0));
    vec4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    vec2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);
    return o4.y * d.y + o4.x * (1.0 - d.y);
}

// Random rotation matrices 
const mat3 rot1 = mat3(-0.37, 0.36, 0.85,-0.14,-0.93, 0.34,0.92, 0.01,0.4);
const mat3 rot2 = mat3(-0.55,-0.39, 0.74, 0.33,-0.91,-0.24,0.77, 0.12,0.63);
const mat3 rot3 = mat3(-0.71, 0.52,-0.47,-0.08,-0.72,-0.68,-0.7,-0.45,0.56);

// Random rotation to reduce banding 
float simplex3d_fractal(vec3 m) 
{
	// weight * noise(m * rotation matrix) 
    return   0.5333333 * noise(m*rot1)
			+ 0.2666667 * noise(2.0*m*rot2)
			+ 0.1333333 * noise(4.0*m*rot3)
			+ 0.0666667 * noise(8.0*m);
}

float ManhattanDist(vec3 a, vec3 b) {
	return abs(a.x - b.x) + abs(a.y - b.y) + abs(a.z - b.z);
}

vec3 GetVolumetricFog(vec3 p, vec3 c) {

	float TotalFog = 0.0f;

	for (int light = 0 ; light < u_LightCount ; light++) {
		
		vec3 LightAt = LightLocations[light].xyz;
		LightAt += vec3(0.5f);
		
		if (ManhattanDist(LightAt, p) < 6.0f) 
		{
			float EuclideanDistance = pow(distance(p,LightAt)*8.0f, 1.0f);

			if (EuclideanDistance < 4.0f * 4.0f) {

				TotalFog += (1.0f / (EuclideanDistance*EuclideanDistance)) * 64.0f;
			}
		}
	}

	return TotalFog * c;
}

void main() 
{
	// Ray properties
	vec3 rO = u_ViewerPosition;

	// Dither
	vec3 BlueNoise = texture(u_BlueNoise, v_TexCoords * (u_Dimensions / vec2(256.0f))).xyz;

	vec3 rD = GetRayDirectionAt(v_TexCoords);
	vec3 TotalLighting = vec3(0.0f);

	float BaseLinearDepth = texture(u_LinearDepthTexture, v_TexCoords).x;
	BaseLinearDepth = BaseLinearDepth < 0.0f ? 10000.0f : BaseLinearDepth;

	vec3 WorldPosition = rO;
	vec3 RayDirection = rD;
	RayDirection = normalize(RayDirection);
	int DensitySamples = 0;

	bool Use3DNoiseForDensity = true;
	
	//Ray march through participating media and gather densities 
	for (int x = 0; x < 60; x++)
	{
		if (!IsInVolume(WorldPosition))
		{
			break;
		}

		// Depth test 
		float DistanceFromCamera = distance(WorldPosition, u_ViewerPosition);
		if (DistanceFromCamera >= BaseLinearDepth) {
			break;
		}
		
		uint Sample = GetVoxel(ivec3(WorldPosition));
		uint Unpacked1 = Sample & 0xFF;
		uint Unpacked2 = (Sample >> 8) & 0xFF;
		int InitialDistance = int(Unpacked1);
		int BlockType = int(Unpacked2);

		if (InitialDistance == 0) {
			WorldPosition += RayDirection * BlueNoise.x;
			continue;
		}

		vec3 Color = BlockAverageColorData[BlockType].xyz;
		TotalLighting += GetVolumetricFog(WorldPosition,vec3(Color));
		WorldPosition += RayDirection * BlueNoise.x;
	}

	TotalLighting *= 2.7525f;
	o_Volumetrics = vec3(TotalLighting / 60.0f);
}