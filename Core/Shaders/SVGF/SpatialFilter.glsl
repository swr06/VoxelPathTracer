#version 330 core

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out float o_Variance;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_SH;
uniform sampler2D u_CoCg;
uniform sampler2D u_Utility;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_BlockIDTexture;
uniform sampler2D u_VarianceTexture;

uniform vec2 u_Dimensions;
uniform int u_Step;
uniform bool u_ShouldDetailWeight;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

const float POSITION_THRESH = 4.0f;

vec3 GetRayDirectionAt(vec2 txc)
{
	vec4 clip = vec4(txc * 2.0f - 1.0f, -1.0f, 1.0f);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0f, 0.0f);
	return vec3(u_InverseView * eye);
}

vec4 GetPositionAt(vec2 txc)
{
	float Dist = texture(u_PositionTexture, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

int GetBlockAt(vec2 txc)
{
	float id = texture(u_BlockIDTexture, txc).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

bool IsSky(float d)
{
	return d < 0.0f;
}

bool InScreenSpace(in vec2 v) 
{
	return v.x > 0.0f && v.x < 1.0f && v.y > 0.0f && v.y < 1.0f;
}

vec3 Saturate(vec3 x)
{
	return clamp(x, 0.0f, 1.0f);
}

float GetVarianceEstimate(out float BaseVariance)
{
	vec2 TexelSize = 1.0f / textureSize(u_SH, 0);
	float VarianceSum = 0.0f;

	const float Kernel[3] = { 1.0 / 4.0, 1.0 / 8.0, 1.0 / 16.0 };

	for (int x = -1 ; x <= 1 ; x++)
	{
		for (int y = -1 ; y <= 1 ; y++)
		{
			vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;
			
			if (!InScreenSpace(SampleCoord)) { continue ; }

			float KernelValue = Kernel[abs(x) + abs(y)]; 
			float V = texture(u_VarianceTexture, SampleCoord).r;

			if (x == 0 && y == 0) { BaseVariance = V; }

			VarianceSum += V * KernelValue;
		}
	}

	return VarianceSum;
}

// From quake2rtx 
vec3 SHToIrridiance(vec4 shY, vec2 CoCg, vec3 v)
{
    float x = dot(shY.xyz, v);
    float Y = 2.0 * (1.023326f * x + 0.886226f * shY.w);
    Y = max(Y, 0.0);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

float sqr(float x) { return x * x; }

void main()
{
	const float AtrousWeights[3] = { 1.0f, 2.0f / 3.0f, 1.0f / 6.0f };
	const float AtrousWeights2[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);

	vec2 TexelSize = 1.0f / u_Dimensions;
	vec4 TotalColor = vec4(0.0f);

	vec4 BasePosition = GetPositionAt(v_TexCoords);
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec3 BaseUtility = texture(u_Utility, v_TexCoords).xyz;
	vec4 BaseSH = texture(u_SH, v_TexCoords).xyzw;
	vec2 BaseCoCg = texture(u_CoCg, v_TexCoords).xy;
	vec3 BaseRadiance = SHToIrridiance(BaseSH, BaseCoCg, BaseNormal);
	float BaseLuminance = GetLuminance(BaseRadiance);
	float BaseVariance = 0.0f;
	float VarianceEstimate = GetVarianceEstimate(BaseVariance);

	// Sigma :

	// Start with the base inputs, one iteration of the loop can then be skipped
	vec4 TotalSH = BaseSH;
	vec2 TotalCoCg = BaseCoCg;
	float TotalWeight = 1.0f;
	float TotalVariance = BaseVariance;

	float PhiColor = 1.0e1 * sqrt(max(0.0f, 1e-10f + VarianceEstimate));

	for (int x = -2 ; x <= 2 ; x++)
	{
		for (int y = -2 ; y <= 2 ; y++)
		{
			vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;

			if (!InScreenSpace(SampleCoord)) { continue; }
			if (x == 0 && y == 0) { continue ; }

			// Sample :
			vec4 SampleSH = texture(u_SH, SampleCoord).xyzw;
			vec2 SampleCoCg = texture(u_CoCg, SampleCoord).xy;
			vec3 SampleNormal = texture(u_NormalTexture, SampleCoord).xyz;
			vec3 SamplePosition = GetPositionAt(SampleCoord).xyz;
			vec3 SampleRadiance = SHToIrridiance(SampleSH, SampleCoCg, SampleNormal);
			float SampleLuma = GetLuminance(SampleRadiance);
			float SampleVariance = texture(u_VarianceTexture, SampleCoord).r;

			// Weights : 
			vec3 PositionDifference = abs(SamplePosition.xyz - BasePosition.xyz);
            float DistSqr = dot(PositionDifference, PositionDifference);

			if (DistSqr < 2.4f) {
				float NormalWeight = pow(max(dot(BaseNormal, SampleNormal), 0.0f), 12.0f);
				float LuminosityWeight = abs(SampleLuma - BaseLuminance) / 1.0e1;
				float Weight = exp(-LuminosityWeight - NormalWeight);
				Weight = max(Weight, 0.0f);

				// Kernel Weights : 
				float XWeight = AtrousWeights[abs(x)];
				float YWeight = AtrousWeights[abs(y)];
				Weight = (XWeight * YWeight) * Weight;

				TotalSH += SampleSH * Weight;
				TotalCoCg += SampleCoCg * Weight;
				TotalVariance += sqr(Weight) * SampleVariance;
				TotalWeight += Weight;
			}
		}
	}
	
	TotalSH /= TotalWeight;
	TotalCoCg /= TotalWeight;
	TotalVariance /= sqr(TotalWeight);
	
	// Output : 
	o_SH = TotalSH;
	o_CoCg = TotalCoCg;
	o_Variance = TotalVariance;
}