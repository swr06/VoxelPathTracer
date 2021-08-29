#version 430 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

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

layout (std430, binding = 2) buffer SSBO_BlockAverageData
{
    vec3 BlockAverageColorData[128];
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

	bool Use3DNoiseForDensity = false;
	
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
		
		// Dither
		vec3 DitheredPosition = WorldPosition.xyz;
		DitheredPosition += BlueNoise * vec3(0.5);
		uint Sample = GetVoxel(ivec3(DitheredPosition));
		uint Unpacked1 = Sample & 0xFF;
		uint Unpacked2 = (Sample >> 8) & 0xFF;
		int InitialDistance = int(Unpacked1);
		int BlockType = int(Unpacked2);

		if (InitialDistance == 0) 
		{ 
			WorldPosition += RayDirection * BlueNoise.x;
			continue; 
		}

		vec3 Color = BlockAverageColorData[BlockType];

		float Lighting = 0.0f; 
		float DistSqr = InitialDistance;
		DistSqr = DistSqr * DistSqr;
		Lighting = DistSqr;

		const float LightingPow = 1.0f / (sqrt(2.0f));
		const float Sqrt2 = sqrt(2.0f);

		Lighting = pow(Lighting, LightingPow);
		
		// Add variation using 3D noise

		if (Use3DNoiseForDensity) {
			float Noise3D = noise(WorldPosition);
			float NoiseFactor = pow(Noise3D, Sqrt2);
			Lighting = Lighting * NoiseFactor; 
		}

		else {

			Lighting *= 0.75750;
		}

		TotalLighting += Color * (Lighting * (1.0f + 0.1f));
		WorldPosition += RayDirection * BlueNoise.x;
	}

	o_Volumetrics = vec3(TotalLighting / 60.0f);
}