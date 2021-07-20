#version 330 core

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out vec3 o_Utility;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousPositionTexture;

uniform sampler2D u_CurrentSH;
uniform sampler2D u_PreviousSH;
uniform sampler2D u_CurrentCoCg;
uniform sampler2D u_PrevCoCg;

uniform sampler2D u_CurrentNormalTexture;
uniform sampler2D u_PreviousNormalTexture;
uniform sampler2D u_NoisyLuminosity;

uniform sampler2D u_PreviousUtility;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_MinimumMix = 0.25f;
uniform float u_MaximumMix = 0.975f;
uniform vec3 u_PrevCameraPos;
uniform vec3 u_CurrentCameraPos;

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

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

bool InScreenSpace(vec2 x)
{
    return x.x < 1.0f && x.x > 0.0f && x.y < 1.0f && x.y > 0.0f;
}

bool InThresholdedScreenSpace(vec2 x)
{
	const float b = 0.01f;
    return x.x < 1.0f - b && x.x > b && x.y < 1.0f - b && x.y > b;
}

void main()
{
	vec4 BasePosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords);
	vec3 BaseNormal = texture(u_CurrentNormalTexture, v_TexCoords).xyz;
	vec4 BaseSH = texture(u_CurrentSH, v_TexCoords).rgba;
	vec2 BaseCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;

	vec2 TexelSize = 1.0f / textureSize(u_CurrentSH, 0);

	float TotalWeight = 0.0f;
	vec4 SumSH = vec4(0.0f); 
	vec2 SumCoCg = vec2(0.0f);
	float SumLuminosity = 0.0f;
	float SumSPP = 0.0f;
	float SumMoment = 0.0f;
	vec2 ReprojectedCoord = Reprojection(BasePosition.xyz);
	vec2 PositionFract = fract(ReprojectedCoord);
	vec2 OneMinusPositionFract = 1.0f - PositionFract;
	float BaseLuminosity = texture(u_NoisyLuminosity, v_TexCoords).r;

	// https://www.eso.org/sci/software/esomidas/doc/user/18NOV/volb/node317.html
	float Weights[5] = float[5](3.0f / 32.0f,
								3.0f / 32.0f,
								9.0f / 64.0f,
								3.0f / 32.0f,
								3.0f / 32.0f);

	const vec2 Offsets[5] = vec2[5](vec2(1, 0), vec2(0, 1), vec2(0.0f), vec2(-1, 0), vec2(0, -1));

	for (int i = 0 ; i < 5 ; i++)
	{
		vec2 Offset = Offsets[i];
		vec2 SampleCoord = ReprojectedCoord + vec2(Offset.x, Offset.y) * TexelSize;

		if (!InThresholdedScreenSpace(SampleCoord)) { continue; }

		vec3 PreviousPositionAt = GetPositionAt(u_PreviousPositionTexture, SampleCoord).xyz;
		vec3 PreviousNormalAt = texture(u_PreviousNormalTexture, SampleCoord).xyz;
		vec3 PositionDifference = abs(BasePosition.xyz - PreviousPositionAt.xyz);
		//float PositionError = dot(PositionDifference, PositionDifference);
		float PositionError = distance(PreviousPositionAt.xyz, BasePosition.xyz);
		float CurrentWeight = Weights[i];

		if (PositionError < 1.5f &&
			PreviousNormalAt == BaseNormal)
		{
			vec3 PreviousUtility = texture(u_PreviousUtility, SampleCoord).xyz;
			vec4 PreviousSH = texture(u_PreviousSH, SampleCoord).xyzw;
			vec2 PreviousCoCg = texture(u_PrevCoCg, SampleCoord).xy;

			SumSH += PreviousSH * CurrentWeight;
			SumCoCg += PreviousCoCg * CurrentWeight;
			SumSPP += PreviousUtility.x * CurrentWeight;
			SumMoment += PreviousUtility.y * CurrentWeight;
			SumLuminosity += PreviousUtility.z * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	if (TotalWeight > 0.0f) { 
		SumSH /= TotalWeight;
		SumCoCg /= TotalWeight;
		SumMoment /= TotalWeight;
		SumSPP /= TotalWeight;
		SumLuminosity /= TotalWeight;
	}

	const bool AccumulateAll = false;

	float BlendFactor = max(1.0f / (SumSPP + 1.0f), 0.05f);
	float MomentFactor = max(1.0f / (SumSPP + 1.0f), 0.05f);

	
	if (AccumulateAll) {
		BlendFactor = 0.01f;
	}

	float UtilitySPP = SumSPP + 1.0;
	float UtilityMoment = (1 - MomentFactor) * SumMoment + MomentFactor * pow(BaseLuminosity, 2.0f);
	
	float CurrentNoisyLuma = texture(u_NoisyLuminosity, v_TexCoords).r;
	float StoreLuma = mix(SumLuminosity, CurrentNoisyLuma, BlendFactor).r;

	o_SH = mix(SumSH, BaseSH, BlendFactor);
	o_CoCg = mix(SumCoCg, BaseCoCg, BlendFactor);
	o_Utility = vec3(UtilitySPP, UtilityMoment, StoreLuma);
}


