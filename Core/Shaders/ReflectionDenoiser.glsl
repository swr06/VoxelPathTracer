#version 330 core

#define THRESH 4.4

layout (location = 0) out vec4 o_SpatialResult;
layout (location = 1) out vec2 o_CoCg;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_InputTexture;
uniform sampler2D u_InputCoCgTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_PBRTex;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform vec2 u_Dimensions;
uniform int u_Step;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevView;
uniform mat4 u_PrevProjection;

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

vec3 ProjectPositionPrevious(vec3 pos)
{
	vec3 WorldPos = pos;
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;

	return ProjectedPosition.xyz;
}

vec2 Reprojection(vec3 pos) 
{
	return ProjectPositionPrevious(pos).xy * 0.5f + 0.5f;
}

// Edge stopping function
bool SampleValid(in vec2 SampleCoord, in vec3 InputPosition, in vec3 InputNormal, in vec3 Col, in float BaseLuminance, in float BaseRoughness, out float SampleDistance)
{
	bool InScreenSpace = SampleCoord.x > 0.0f && SampleCoord.x < 1.0f && SampleCoord.y > 0.0f && SampleCoord.y < 1.0f;
	vec4 PositionAt = GetPositionAt(u_PositionTexture, SampleCoord);
	vec3 NormalAt = texture(u_NormalTexture, SampleCoord).xyz;
	SampleDistance = distance(PositionAt.xyz, InputPosition.xyz);
	SampleDistance = clamp(SampleDistance, 0.0f, THRESH);
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
	return vec4(clamp(Col.xyz, 0.0f, 1.0f + 0.8f), Col.w);
}

void main()
{
	vec4 BlurredColor = vec4(0.0f);
	vec2 BlurredCoCg = vec2(0.0f);
	vec3 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec4 BaseColor = texture(u_InputTexture, v_TexCoords).xyzw;
	float BaseLuminance = GetLuminance(BaseColor.xyz);

	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	vec2 Reprojected = Reprojection(BasePosition);
	float RoughnessAt;
	float ReprojectionBias = 0.01f;
	
	if (Reprojected.x > 0.0f + ReprojectionBias && Reprojected.x < 1.0f - ReprojectionBias &&
		Reprojected.y > 0.0f + ReprojectionBias && Reprojected.y < 1.0f - ReprojectionBias) {
		RoughnessAt = texture(u_PBRTex, Reprojected).r;
	} else { RoughnessAt = 0.01f; }

	for (int s = 0 ; s < GAUSS_KERNEL; s++)
	{
		float CurrentWeight = GaussianWeightsNormalized[s];
		int Sample = GaussianOffsets[s]; // todo : use u_Step here!
		
		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (Sample * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (Sample * TexelSize));
		vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;
		vec2 SampleCoCg = texture(u_InputCoCgTexture, SampleCoord).rg;

		// Luminosity weights
		float LumaAt = GetLuminance(SampleColor.xyz);
		float LumaDifference = abs(BaseLuminance - LumaAt);
		float LumaThreshold = RoughnessAt * (1.0f + RoughnessAt);
		LumaThreshold = clamp(LumaThreshold * 0.750, 0.01f, 100.0f);

		float SampleDistance = 0.0f;
		
		if (SampleValid(SampleCoord, BasePosition, BaseNormal, SampleColor.xyz, BaseLuminance, RoughnessAt, SampleDistance) && LumaDifference < LumaThreshold)
		{
			float DistanceWeight = abs(THRESH - SampleDistance);
			DistanceWeight = pow(DistanceWeight, 5.0f);
			//vec3 SampledNormalMapped = texture(u_NormalMappedTexture, SampleCoord).xyz;
			//float NormalWeight = pow(abs(dot(SampledNormalMapped, BaseNormalMapped)), 16.0f);
			float NormalWeight = 1.0f;

			CurrentWeight *= NormalWeight * DistanceWeight;
			BlurredColor += SampleColor * CurrentWeight;
			BlurredCoCg += SampleCoCg * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	BlurredColor = BlurredColor / max(TotalWeight, 0.01f);
	BlurredCoCg = BlurredCoCg / max(TotalWeight, 0.01f);

	vec4 BaseSampled = texture(u_InputTexture, v_TexCoords).rgba;
	vec2 BaseSampledCoCg = texture(u_InputCoCgTexture, v_TexCoords).rg;

	float m = 1.0f;

	// Preserve a few more details in smoother materials
	//if (RoughnessAt <= 0.125f)
	//{
	//	m = RoughnessAt * 16.0f;
	//	m = 1.0f - m;
	//	m = pow(m, 6.0f);
	//	m = clamp(m, 0.01f, 0.999f);
	//} 

	if (RoughnessAt <= 0.125f)
	{
		m = 0.125f - RoughnessAt;
		m = pow(m, 8.0f);
		m = clamp(m, 0.01f, 0.999f);
	} 


	o_SpatialResult = mix(BaseSampled, BlurredColor, m);
	o_CoCg = mix(BaseSampledCoCg, BlurredCoCg, m);
}