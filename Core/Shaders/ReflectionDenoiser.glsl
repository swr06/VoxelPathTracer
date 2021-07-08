#version 330 core

#define THRESH 0.3f

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_NormalMappedTexture;
uniform sampler2D u_PBRTex;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform vec2 u_Dimensions;
uniform int u_Step;

// Large kernel gaussian denoiser //
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

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(v_RayDirection) * Dist, Dist);
}

// Edge stopping function
bool SampleValid(in vec2 SampleCoord, in vec3 InputPosition, in vec3 InputNormal, in vec3 Col, in float BaseLuminance, in float BaseRoughness)
{
	bool InScreenSpace = SampleCoord.x > 0.0f && SampleCoord.x < 1.0f && SampleCoord.y > 0.0f && SampleCoord.y < 1.0f;
	vec4 PositionAt = GetPositionAt(u_PositionTexture, SampleCoord);
	vec3 NormalAt = texture(u_NormalTexture, SampleCoord).xyz;
	return (abs(PositionAt.z - InputPosition.z) <= THRESH) 
			&& (abs(PositionAt.x - InputPosition.x) <= THRESH) 
			&& (abs(PositionAt.y - InputPosition.y) <= THRESH) 
			&& (InputNormal == NormalAt) 
			&& (PositionAt.w > 0.0f) && (InScreenSpace);
			//&& pow(abs(BaseLuminance - Luma), 1.0f) < (BaseRoughness);
}



// Basic clamp firefly reject
vec4 FireflyReject(vec4 Col)
{
	return vec4(clamp(Col.xyz, 0.0f, 1.0f + 1.666f), Col.w);
}

void main()
{
	vec4 BlurredColor = vec4(0.0f);
	vec3 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec3 BaseColor = texture(u_InputTexture, v_TexCoords).xyz;
	float BaseLuminance = GetLuminance(BaseColor);

	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	float RoughnessAt = texture(u_PBRTex, v_TexCoords).r;
	vec3 BaseNormalMapped = texture(u_NormalMappedTexture, v_TexCoords).rgb;

	for (int s = 0 ; s < GAUSS_KERNEL; s++)
	{
		float CurrentWeight = GaussianWeightsNormalized[s];
		int Sample = GaussianOffsets[s]; // todo : use u_Step here!
		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (Sample * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (Sample * TexelSize));
		vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;
		float LumaAt = GetLuminance(SampleColor.xyz);
		float LumaDifference = abs(BaseLuminance - LumaAt);
		float LumaThreshold = RoughnessAt * (1.0f + RoughnessAt);
		LumaThreshold = clamp(LumaThreshold * 0.750, 0.01f, 100.0f);
		//LumaThreshold = exp(RoughnessAt + 1);
		
		if (SampleValid(SampleCoord, BasePosition, BaseNormal, SampleColor.xyz, BaseLuminance, RoughnessAt) && LumaDifference < LumaThreshold)
		{
			vec3 SampledNormalMapped = texture(u_NormalMappedTexture, SampleCoord).xyz;
			float NormalWeight = pow(abs(dot(SampledNormalMapped, BaseNormalMapped)), 16.0f);
			CurrentWeight *= NormalWeight;
			BlurredColor += FireflyReject(SampleColor) * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	BlurredColor = BlurredColor / max(TotalWeight, 0.01f);

	vec4 BaseSampled = texture(u_InputTexture, v_TexCoords).rgba;
	float m = 1.0f;

	// Preserve a few more details in smoother materials
	if (RoughnessAt <= 0.125f)
	{
		m = RoughnessAt * 16.0f;
		m = 1.0f - m;
		m = clamp(m, 0.01f, 0.999f);
	}

	o_SpatialResult = mix(BaseSampled, BlurredColor, m);
}