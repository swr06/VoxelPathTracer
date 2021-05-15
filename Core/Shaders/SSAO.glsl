#version 330 core

layout (location = 0) out float o_AOValue;

in vec2 v_TexCoords;

uniform sampler2D u_DepthTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_NoiseTexture;
uniform sampler2D u_SSAOKernel;
uniform vec2 u_Dimensions;

uniform mat4 u_InverseViewMatrix;
uniform mat4 u_InverseProjectionMatrix;
uniform mat4 u_ViewMatrix;
uniform mat4 u_ProjectionMatrix;

uniform int SAMPLE_SIZE = 28;

float CalcViewZ(vec2 Coords)
{
    float Depth = texture(u_DepthTexture, Coords).x;
    float ViewZ = u_ProjectionMatrix[3][2] / (2 * Depth -1  - u_ProjectionMatrix[2][2]);
    return ViewZ;
}

vec2 NoiseScale = vec2(u_Dimensions.x / 8.0f, u_Dimensions.y / 8.0f);

vec2 Convert1DTo2DIndex(in float idx, in vec2 texsize)
{
	return vec2(
		idx / texsize.x,
		mod(idx, texsize.y)
	);
}

vec3 WorldPosFromDepth(float depth) 
{
    float z = depth * 2.0 - 1.0;

    vec4 clipSpacePosition = vec4(v_TexCoords * 2.0 - 1.0, z, 1.0);
    vec4 viewSpacePosition = u_InverseProjectionMatrix * clipSpacePosition;

    // Perspective division
    viewSpacePosition /= viewSpacePosition.w;
    vec4 worldSpacePosition = u_InverseViewMatrix * viewSpacePosition;
    return worldSpacePosition.xyz;
}

vec3 ViewPosFromDepth(float depth) 
{
    float z = depth * 2.0 - 1.0;

    vec4 clipSpacePosition = vec4(v_TexCoords * 2.0 - 1.0, z, 1.0);
    vec4 viewSpacePosition = u_InverseProjectionMatrix * clipSpacePosition;
    viewSpacePosition /= viewSpacePosition.w;

    return viewSpacePosition.xyz;
}

float CompareDepthsFar(float depth1, float depth2) 
{
	float garea = 2.0; 
	float diff = (depth1 - depth2) * 100.0; 

	if (diff < 0.4f)
	{
		return 0.0f;
	} 
	
	else 
	{
		return 1.0f;
	}
}

const float Radius = 0.69f; //nice
const float Bias = 0.14f;

void main()
{
	float PixelDepth = texture(u_DepthTexture, v_TexCoords).r;

	if (PixelDepth == 1.0f)
	{
		o_AOValue = 0.0f;
		return;
	}

	vec3 Position = ViewPosFromDepth(PixelDepth);

	vec3 NoiseValue = normalize(vec3(texture(u_NoiseTexture, v_TexCoords * NoiseScale).rg, 0.0f));
	vec3 Normal = normalize(vec3(u_ViewMatrix * vec4(texture(u_NormalTexture, v_TexCoords).xyz, 0.0f)));
	vec3 Tangent = normalize(NoiseValue - Normal * dot(NoiseValue, Normal));
	vec3 Bitangent = normalize(cross(Normal, Tangent));
	mat3 TBN = mat3(Tangent, Bitangent, Normal);

	vec2 KernelTextureSize = textureSize(u_SSAOKernel, 0);
	int KernelSampleCount = int(KernelTextureSize.x) * int(KernelTextureSize.y);
	o_AOValue = 0.0f;


	for (int i = 0 ; i < clamp(SAMPLE_SIZE, 1, KernelSampleCount) ; i++)
	{
		vec3 KernelValue = texelFetch(u_SSAOKernel, ivec2(Convert1DTo2DIndex(i, KernelTextureSize)), 0).rgb;
		vec3 OrientedKernel = TBN * KernelValue;

		vec3 SamplePosition = Position + (OrientedKernel * Radius);
		vec4 ProjectedPosition = u_ProjectionMatrix * vec4(SamplePosition, 1.0f);
		ProjectedPosition.xyz /= ProjectedPosition.w; // Perspective division
		ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;

		if (ProjectedPosition.x > 0.0f && ProjectedPosition.y > 0.0f && ProjectedPosition.x < 1.0f && ProjectedPosition.y < 1.0f)
		{
			float NonLinearSampleDepth = texture(u_DepthTexture, ProjectedPosition.xy).r;
			float SampleDepth = ViewPosFromDepth(NonLinearSampleDepth).z;

			float RangeCheck = smoothstep(0.0, 1.0, Radius / abs(Position.z - SampleDepth));

			//o_AOValue += (SampleDepth >= SamplePosition.z + Bias ? 1.0 : 0.0) * RangeCheck; 
			o_AOValue += (CompareDepthsFar(SampleDepth, SamplePosition.z + Bias)) * RangeCheck;
		}
	}

	o_AOValue /= SAMPLE_SIZE;
	o_AOValue = 1.0f - o_AOValue;
}

