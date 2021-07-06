#version 330 core

#define THRESH 1.41421f

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;

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

// Edge stopping function
bool SampleValid(in vec2 SampleCoord, in vec3 InputPosition, in vec3 InputNormal)
{
	bool InScreenSpace = SampleCoord.x > 0.0f && SampleCoord.x < 1.0f && SampleCoord.y > 0.0f && SampleCoord.y < 1.0f;
	vec4 PositionAt = texture(u_PositionTexture, SampleCoord);
	vec3 NormalAt = texture(u_NormalTexture, SampleCoord).xyz;
	return (abs(PositionAt.z - InputPosition.z) <= THRESH) 
			&& (abs(PositionAt.x - InputPosition.x) <= THRESH) 
			&& (abs(PositionAt.y - InputPosition.y) <= THRESH) 
			&& (InputNormal == NormalAt) 
			&& (PositionAt.w > 0.0f) && (InScreenSpace);
}

void main()
{
	vec4 BlurredColor = vec4(0.0f);
	vec3 BasePosition = texture(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;

	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	int SZ = u_LargeKernel ? GAUSS_KERNEL : GAUSS_KERNEL_SMALL;

	for (int s = 0 ; s < SZ ; s++)
	{
		float CurrentWeight = u_LargeKernel ? GaussianWeightsNormalized[s] : GaussianWeightsNormalized_SMALL[s];
		int Sample = u_LargeKernel ? GaussianOffsets[s] : GaussianOffsets_SMALL[s]; 

		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (Sample * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (Sample * TexelSize));

		if (SampleValid(SampleCoord, BasePosition, BaseNormal))
		{
			vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;
			BlurredColor += SampleColor * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	BlurredColor = BlurredColor / max(TotalWeight, 0.01f);

	o_SpatialResult = BlurredColor;
}