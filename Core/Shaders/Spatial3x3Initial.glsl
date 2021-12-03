#version 330 core

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out float o_Utility;
layout (location = 3) out vec2 o_AOSky;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_SH;
uniform sampler2D u_CoCg;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_AO;
uniform sampler2D u_Utility;

uniform vec2 u_Dimensions;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

float u_ColorPhiBias = 2.66f;
uniform float u_DeltaTime;
uniform float u_Time;

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

float SHToY(vec4 shY)
{
    return max(0, 3.544905f * shY.w);
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

vec3 SampleNormalFromTex(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}

float sqr(float x) { return x * x; }
float GetSaturation(in vec3 v) { return length(v); }

float GradientNoise()
{
	vec2 coord = gl_FragCoord.xy + mod(u_Time * 100.493850275f, 500.0f);
	float noise = fract(52.9829189f * fract(0.06711056f * coord.x + 0.00583715f * coord.y));
	return noise;
}

void main()
{
	const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );
	const float AtrousWeights2[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);

	vec2 TexelSize = 1.0f / u_Dimensions;
	vec4 TotalColor = vec4(0.0f);


	vec4 BasePosition = GetPositionAt(v_TexCoords);
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;
	vec4 BaseSH = texture(u_SH, v_TexCoords).xyzw;
	vec2 BaseCoCg = texture(u_CoCg, v_TexCoords).xy;
	float BaseLuminance = SHToY(BaseSH);
	vec2 BaseAOSky = texture(u_AO, v_TexCoords).xy;
	float BaseUtility = texture(u_Utility, v_TexCoords).x;

	// Start with the base inputs, one iteration of the loop can then be skipped
	vec4 TotalSH = BaseSH;
	vec2 TotalCoCg = BaseCoCg;
	//float LumaTotal = BaseUtility;
	float TotalWeight = 1.0f;
	vec2 TotalAOSky = BaseAOSky;
	float TotalAOWeight = 1.0f;
	float PhiColor = 4.0f;
	ivec2 Jitter = ivec2((GradientNoise() - 0.5f) * float(1.0f * 1.1f));

	for (int x = -1 ; x <= 1 ; x++)
	{
		for (int y = -1 ; y <= 1 ; y++)
		{
			if (x == 0 && y == 0) { continue ; }
			vec2 SampleCoord = v_TexCoords + ((vec2(x, y) * float(1))) * TexelSize;
			if (!InScreenSpace(SampleCoord)) { continue; }

			vec3 SamplePosition = GetPositionAt(SampleCoord).xyz;

			// Weights : 
			vec3 PositionDifference = abs(SamplePosition.xyz - BasePosition.xyz);
            float DistSqr = dot(PositionDifference, PositionDifference);

			if (DistSqr < 1.0f) {

				vec4 SampleSH = texture(u_SH, SampleCoord).xyzw;
				vec2 SampleCoCg = texture(u_CoCg, SampleCoord).xy;
				vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
				float SampleLuma = SHToY(SampleSH);
				float NormalWeight = pow(max(dot(BaseNormal, SampleNormal), 0.0f), 16.0f);
				float LuminosityWeight = abs(SampleLuma - BaseLuminance) / PhiColor;
				float Weight = exp(-LuminosityWeight - NormalWeight);
				Weight = max(Weight, 0.01f);
				float XWeight = AtrousWeights[abs(x)];
				float YWeight = AtrousWeights[abs(y)];
				Weight = (XWeight * YWeight) * Weight;
				Weight = max(Weight, 0.01f);
				Weight = clamp(Weight, 0.0f, 1.0f);
				
				//LumaTotal += texture(u_Utility,SampleCoord).x * Weight;
				TotalSH += SampleSH * Weight;
				TotalCoCg += SampleCoCg * Weight;
				TotalWeight += Weight;
				TotalAOSky += texture(u_AO, SampleCoord).xy * Weight;
				TotalAOWeight += Weight;
			}
		}
	}
	
	TotalWeight = max(TotalWeight, 0.01f);
	TotalSH /= TotalWeight;
	TotalCoCg /= TotalWeight;
	TotalAOSky /= max(TotalAOWeight, 0.01f);
	//LumaTotal /= TotalWeight;
	
	// Output : 
	o_SH = TotalSH;
	o_CoCg = TotalCoCg;
	o_AOSky = TotalAOSky;

	//o_Utility = LumaTotal;
	o_Utility = BaseUtility;
}