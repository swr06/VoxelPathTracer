#version 330 core

#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out vec3 o_Utility;
layout (location = 3) out vec2 o_AOAndSkyLighting;

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

uniform sampler2D u_CurrentBlockIDTexture;
uniform sampler2D u_PrevBlockIDTexture;

uniform sampler2D u_CurrentAO;
uniform sampler2D u_PreviousAO;

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
uniform bool u_BeUseful;

uniform float u_DeltaTime;
uniform float u_Time;

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
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
	const float b = 0.0035f;
    return x.x < 1.0f - b && x.x > b && x.y < 1.0f - b && x.y > b;
}

bool CompareFloatNormal(float x, float y) {
    return abs(x - y) < 0.02f;
}

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

vec3 SampleNormal(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}

float GradientNoise()
{
	vec2 coord = gl_FragCoord.xy + mod(u_Time * 100.493850275f, 500.0f);
	float noise = fract(52.9829189f * fract(0.06711056f * coord.x + 0.00583715f * coord.y));
	return noise;
}

void main()
{
	vec4 BasePosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords);
	vec3 BaseNormal = SampleNormal(u_CurrentNormalTexture, v_TexCoords).xyz;
	vec4 BaseSH = texture(u_CurrentSH, v_TexCoords).rgba;
	vec2 BaseCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;
	vec2 BaseAOSky = texture(u_CurrentAO, v_TexCoords).xy;

	vec2 TexelSize = 1.0f / textureSize(u_CurrentSH, 0);

	float TotalWeight = 0.0f;
	vec4 SumSH = vec4(0.0f); 
	vec2 SumCoCg = vec2(0.0f);
	vec2 SumAOSky = vec2(0.0f);
	float SumLuminosity = 0.0f;
	float SumSPP = 0.0f;
	float SumMoment = 0.0f;
	vec2 ReprojectedCoord = Reprojection(BasePosition.xyz);
	vec2 PositionFract = fract(ReprojectedCoord);
	vec2 OneMinusPositionFract = 1.0f - PositionFract;
	float BaseLuminosity = texture(u_NoisyLuminosity, v_TexCoords).r;

	// Atrous weights : 
	// https://www.eso.org/sci/software/esomidas/doc/user/18NOV/volb/node317.html
	float Weights[5] = float[5](3.0f / 32.0f,
								3.0f / 32.0f,
								9.0f / 64.0f,
								3.0f / 32.0f,
								3.0f / 32.0f);

	const vec2 Offsets[5] = vec2[5](vec2(1, 0), vec2(0, 1), vec2(0.0f), vec2(-1, 0), vec2(0, -1));

	int BaseBlock = clamp(int(floor((texture(u_CurrentBlockIDTexture, v_TexCoords).r) * 255.0f)), 0, 127);
	ivec2 Jitter = ivec2((GradientNoise() - 0.5f) * float(1.5f));

	// Sample neighbours and hope to find a good sample : 
	// Verified that this is indeed *much* better than just sampling at the reprojected coordinate 

	for (int i = 0 ; i < 5 ; i++)
	{
		vec2 Offset = Offsets[i];
		vec2 SampleCoord = ReprojectedCoord + (vec2(Offset.x, Offset.y) + vec2(Jitter)) * TexelSize;
		if (!InThresholdedScreenSpace(SampleCoord)) { continue; }

		vec3 PreviousPositionAt = GetPositionAt(u_PreviousPositionTexture, SampleCoord).xyz;
		vec3 PreviousNormalAt = SampleNormal(u_PreviousNormalTexture, SampleCoord).xyz;
		vec3 PositionDifference = abs(BasePosition.xyz - PreviousPositionAt.xyz);
		float PositionError = dot(PositionDifference, PositionDifference);
		float CurrentWeight = Weights[i];
		float idat = texture(u_PrevBlockIDTexture, SampleCoord).r;
		int SampleBlock = clamp(int(floor((idat) * 255.0f)), 0, 127);

		if (PositionError < 0.8f &&
			PreviousNormalAt == BaseNormal &&
			BaseBlock == SampleBlock)
		{
			vec3 PreviousUtility = texture(u_PreviousUtility, SampleCoord).xyz;
			vec4 PreviousSH = texture(u_PreviousSH, SampleCoord).xyzw;
			vec2 PreviousCoCg = texture(u_PrevCoCg, SampleCoord).xy;

			SumSH += PreviousSH * CurrentWeight;
			SumCoCg += PreviousCoCg * CurrentWeight;
			SumSPP += PreviousUtility.x * CurrentWeight;
			SumMoment += PreviousUtility.y * CurrentWeight;
			SumLuminosity += PreviousUtility.z * CurrentWeight;
			SumAOSky += texture(u_PreviousAO, SampleCoord).xy * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	if (TotalWeight > 0.0f) { 
		SumSH /= TotalWeight;
		SumCoCg /= TotalWeight;
		SumMoment /= TotalWeight;
		SumSPP /= TotalWeight;
		SumLuminosity /= TotalWeight;
		SumAOSky /= TotalWeight;
	}

	const bool AccumulateAll = false;

	float Increment = float(u_BeUseful);
	float BlendFactor = max(1.0f / (SumSPP + Increment), 0.05f);
	float MomentFactor = max(1.0f / (SumSPP + Increment), 0.05f);

	
	if (!u_BeUseful) {
		BlendFactor = 0.99f;
		SumSPP = 0.0f;
	}

	SumSPP = clamp(SumSPP, 0.0f, 120.0f);
	float UtilitySPP = SumSPP + Increment;
	float UtilityMoment = (1.0f - MomentFactor) * SumMoment + MomentFactor * (BaseLuminosity * BaseLuminosity) ;//pow(BaseLuminosity, 2.0f); 
	
	float CurrentNoisyLuma = texture(u_NoisyLuminosity, v_TexCoords).r;
	float StoreLuma = mix(SumLuminosity, CurrentNoisyLuma, BlendFactor);
	
	o_SH = mix(SumSH, BaseSH, BlendFactor);
	o_CoCg = mix(SumCoCg, BaseCoCg, BlendFactor);
	o_AOAndSkyLighting = mix(SumAOSky, BaseAOSky, BlendFactor);
	o_Utility = vec3(UtilitySPP, UtilityMoment, StoreLuma);
	
	// clamp
	o_SH = clamp(o_SH, -100.0f, 100.0f);
	o_CoCg = clamp(o_CoCg, -10.0f, 100.0f);
	o_AOAndSkyLighting = clamp(o_AOAndSkyLighting, 0.0f, 1.0f);
	o_Utility = clamp(o_Utility, -150.0f, 150.0f);
}


