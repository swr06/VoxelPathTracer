#version 330 core
#define PI 3.14159265359

layout (location = 0) out float o_AOValue;

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform vec2 u_Dimensions;

uniform mat4 u_ViewMatrix;
uniform mat4 u_ProjectionMatrix;

uniform float u_Time;

uint SAMPLE_SIZE = 16u;

vec3 ToViewSpace(in vec3 WorldPosition)
{
	return (u_ViewMatrix * vec4(WorldPosition, 1.0f)).xyz;
}

vec2 hammersley2d(uint idx, uint num) 
{
	uint bits = idx;
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	float radicalInverse_VdC = float(bits) * 2.3283064365386963e-10; // / 0x100000000

	return vec2(float(idx) / float(num), radicalInverse_VdC);
}

vec3 hemispherepoint_uniform(float u, float v)
{
	float phi = v * 2 * PI;
	float cosTheta = 1 - u;
	float sinTheta = sqrt(1 - cosTheta * cosTheta);
	return vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

float nextFloat(inout int seed, in float min, in float max);

const float Radius = 0.5f; 
const float Bias = 0.1f;

int RNG_SEED;

void main()
{
	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x);

	vec4 InitialTracePosition = texture(u_PositionTexture, v_TexCoords).rgba;

	if (InitialTracePosition.a <= 0.0f)
	{
		o_AOValue = 0.0f;
		return;
	}

	vec3 Position = ToViewSpace(InitialTracePosition.xyz);
	vec3 Normal = normalize(vec3(u_ViewMatrix * vec4(texture(u_NormalTexture, v_TexCoords).xyz, 0.0f)));
	
	vec3 noise = vec3(nextFloat(RNG_SEED, -1.0f, 1.0f), nextFloat(RNG_SEED, -1.0f, 1.0f), nextFloat(RNG_SEED, -1.0f, 1.0f));
	vec3 Tangent = normalize(noise - Normal * dot(noise, Normal));
	vec3 Bitangent = cross(Normal, Tangent);
	mat3 TBN = mat3(
		normalize(Tangent),
		normalize(Bitangent),
		normalize(Normal));

	o_AOValue = 0.0f;

	for (uint i = 0u ; i < clamp(SAMPLE_SIZE, 1u, 64u) ; i++)
	{
		vec2 Hammersley = hammersley2d(i, SAMPLE_SIZE);
		vec3 Hemisphere = hemispherepoint_uniform(Hammersley.x, Hammersley.y);
		vec3 OrientedKernel = TBN * Hemisphere;

		vec3 SamplePosition = Position + (OrientedKernel * Radius);
		vec4 ProjectedPosition = u_ProjectionMatrix * vec4(SamplePosition, 1.0f);
		ProjectedPosition.xyz /= ProjectedPosition.w; // Perspective division
		ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;

		if (ProjectedPosition.x > 0.0f && ProjectedPosition.y > 0.0f && ProjectedPosition.x < 1.0f && ProjectedPosition.y < 1.0f)
		{
			vec4 SampledPosition = texture(u_PositionTexture, ProjectedPosition.xy).xyzw;
			
			if (SampledPosition.w > 0.0f)
			{
				float SampleDepth = ToViewSpace(SampledPosition.xyz).z;
				float RangeFix = 1.0f - clamp(abs(Position.z - SampleDepth) * 0.6f, 0.0f, 1.0f);

				o_AOValue += (SampleDepth >= SamplePosition.z + Bias ? 1.0 : 0.0)  * RangeFix; 
			}
		}
	}

	o_AOValue /= SAMPLE_SIZE;
	o_AOValue = 1.0f - o_AOValue;
	o_AOValue = clamp(o_AOValue, 0.0f, 1.0f);
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