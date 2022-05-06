#version 330 core

#define PI 3.14159265359
#define TAU (3.14159265359 * 2.0)

layout (location = 0) out float o_AOValue;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_BlueNoise;
uniform vec2 u_Dimensions;

uniform mat4 u_ViewMatrix;
uniform mat4 u_ProjectionMatrix;

uniform float u_SSAOStrength;

uniform float u_Time;
uniform int u_Frame;



uint SAMPLE_SIZE = 24u;


vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(v_RayDirection) * Dist, Dist);
}

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

const float RadiusR = 0.8f; 
const float Bias = 0.21f;

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

vec3 SampleNormalFromTex(sampler2D samp, vec2 txc) { 
    return GetNormalFromID(texture(samp, txc).x);
}

void main()
{
	vec4 InitialTracePosition = GetPositionAt(u_PositionTexture, v_TexCoords).rgba;

	if (InitialTracePosition.a <= 0.0f)
	{
		o_AOValue = 0.0f;
		return;
	}

	vec3 Position = ToViewSpace(InitialTracePosition.xyz);
	vec3 Normal = normalize(vec3(u_ViewMatrix * vec4(SampleNormalFromTex(u_NormalTexture, v_TexCoords), 0.0f)));
	
	vec3 Hash;
	int n = u_Frame % 1024;
	vec2 off = fract(vec2(n * 12664745, n * 9560333) / 16777216.0) * 1024.0;
	ivec2 TextureSize = textureSize(u_BlueNoise, 0);
	ivec2 SampleTexelLoc = ivec2(gl_FragCoord.xy + ivec2(floor(off))) % TextureSize;
	Hash = texelFetch(u_BlueNoise, SampleTexelLoc, 0).xyz;

	vec3 Tangent = normalize(Hash - Normal * dot(Hash, Normal));
	vec3 Bitangent = cross(Normal, Tangent);
	mat3 TBN = mat3(
		normalize(Tangent),
		normalize(Bitangent),
		normalize(Normal));

	o_AOValue = 0.0f;

	float Radius = RadiusR * Hash.y;

	for (uint i = 0u ; i < clamp(SAMPLE_SIZE, 1u, 64u) ; i++)
	{
		vec2 Hammersley = hammersley2d(i, SAMPLE_SIZE);
		vec3 Hemisphere = hemispherepoint_uniform(Hammersley.x, Hammersley.y);
		vec3 OrientedKernel = TBN * Hemisphere;

		vec3 SamplePosition = Position + (OrientedKernel * Radius);
		vec4 ProjectedPosition = u_ProjectionMatrix * vec4(SamplePosition, 1.0f);
		ProjectedPosition.xyz /= ProjectedPosition.w; // Perspective division
		ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;

		if (ProjectedPosition.xy == clamp(ProjectedPosition.xy, 0.0001, 0.9999f))
		{
			vec4 SampledPosition = GetPositionAt(u_PositionTexture, ProjectedPosition.xy).xyzw;
			
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
	o_AOValue = clamp(pow(o_AOValue, clamp(u_SSAOStrength, 0.01f, 4.0f)), 0.0f, 1.0f);
}

