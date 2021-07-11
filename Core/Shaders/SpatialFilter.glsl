#version 330 core

#define THRESH 4.0f

/*
-- KNOWN ISSUES --
- Overblurs TOO much
- Uses a gaussian kernel, which isnt ideal, atrous would be better and faster 
*/

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_NormalMappedTexture;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform vec2 u_Dimensions;
uniform int u_Step;
uniform bool u_LargeKernel = false;

const int GAUSS_KERNEL = 33;
const float GaussianWeightsNormalized[GAUSS_KERNEL] = float[GAUSS_KERNEL](
	0.004013,
	0.005554,
	0.007527,
	0.00999,
	0.012984,
	0.016524,
	0.020594,
	0.025133,
	0.030036,
	0.035151,
	0.040283,
	0.045207,
	0.049681,
	0.053463,
	0.056341,
	0.058141,
	0.058754,
	0.058141,
	0.056341,
	0.053463,
	0.049681,
	0.045207,
	0.040283,
	0.035151,
	0.030036,
	0.025133,
	0.020594,
	0.016524,
	0.012984,
	0.00999,
	0.007527,
	0.005554,
	0.004013
);

const int GaussianOffsets[GAUSS_KERNEL] = int[GAUSS_KERNEL](
	-16,
	-15,
	-14,
	-13,
	-12,
	-11,
	-10,
	-9,
	-8,
	-7,
	-6,
	-5,
	-4,
	-3,
	-2,
	-1,
	0,
	1,
	2,
	3,
	4,
	5,
	6,
	7,
	8,
	9,
	10,
	11,
	12,
	13,
	14,
	15,
	16
);

const int GAUSS_KERNEL_SMALL = 17;
const float GaussianWeightsNormalized_SMALL[GAUSS_KERNEL_SMALL] = float[GAUSS_KERNEL_SMALL](
	0.030036,
	0.035151,
	0.040283,
	0.045207,
	0.049681,
	0.053463,
	0.056341,
	0.058141,
	0.058754,
	0.058141,
	0.056341,
	0.053463,
	0.049681,
	0.045207,
	0.040283,
	0.035151,
	0.030036
);

const int GaussianOffsets_SMALL[GAUSS_KERNEL_SMALL] = int[GAUSS_KERNEL_SMALL](
	-8,
	-7,
	-6,
	-5,
	-4,
	-3,
	-2,
	-1,
	0,
	1,
	2,
	3,
	4,
	5,
	6,
	7,
	8
);

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

// Edge stopping function
bool SampleValid(in vec2 SampleCoord, in vec3 InputPosition, in vec3 InputNormal, out float Distance)
{
	bool InScreenSpace = SampleCoord.x > 0.0f && SampleCoord.x < 1.0f && SampleCoord.y > 0.0f && SampleCoord.y < 1.0f;
	vec4 PositionAt = GetPositionAt(u_PositionTexture, SampleCoord);
	vec3 NormalAt = texture(u_NormalTexture, SampleCoord).xyz;
	Distance = distance(InputPosition.xyz, PositionAt.xyz);
	Distance = clamp(Distance, 0.0f, THRESH);
	return (abs(PositionAt.z - InputPosition.z) <= THRESH) 
			&& (abs(PositionAt.x - InputPosition.x) <= THRESH) 
			&& (abs(PositionAt.y - InputPosition.y) <= THRESH) 
			&& (InputNormal == NormalAt) 
			&& (PositionAt.w > 0.0f) && (InScreenSpace);
}

vec3 Saturate(vec3 x)
{
	return clamp(x, 0.0f, 1.0f);
}

void main()
{
	vec4 BlurredColor = vec4(0.0f);
	vec3 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec3 BaseColor = texture(u_InputTexture, v_TexCoords).rgb;
	vec3  BaseSampleNormalized = normalize(BaseColor);
    float BaseSaturation  = length(Saturate(BaseColor));

	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	int SZ = u_LargeKernel ? GAUSS_KERNEL : GAUSS_KERNEL_SMALL;

	const float SATURATION_POW = 8.0f;
	const float DOT_POW = 24.0f;

	for (int s = 0 ; s < SZ ; s++)
	{
		float CurrentWeight = u_LargeKernel ? GaussianWeightsNormalized[s] : GaussianWeightsNormalized_SMALL[s];
		int Sample = u_LargeKernel ? GaussianOffsets[s] : GaussianOffsets_SMALL[s]; 

		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (Sample * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (Sample * TexelSize));
		float SamplePositionWeight = 0.0f;
		vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;

		if (SampleValid(SampleCoord, BasePosition, BaseNormal, SamplePositionWeight))
		{
			float DistanceWeight = abs(THRESH - SamplePositionWeight);
			DistanceWeight = pow(DistanceWeight, 6.0f);
			CurrentWeight *= DistanceWeight;
			
			vec3 CurrentSampleNormalized = normalize(SampleColor.xyz);
			float NormalizeWeight = pow(0.5f + 0.5f * dot(BaseSampleNormalized, CurrentSampleNormalized), DOT_POW);
			float SaturationWeight = pow(1.0f - abs(length(Saturate(SampleColor.xyz)) - BaseSaturation), SATURATION_POW);
			//CurrentWeight *= NormalizeWeight * SaturationWeight;
			
			BlurredColor += SampleColor * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	BlurredColor = BlurredColor / max(TotalWeight, 0.01f);
	o_SpatialResult = BlurredColor;
	//o_SpatialResult = texture(u_InputTexture, v_TexCoords).rgba;
}