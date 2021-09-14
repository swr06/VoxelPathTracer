#version 430 core

#define ESTIMATE_VARIANCE_BASED_ON_NEIGHBOURS

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out float o_Variance;
layout (location = 3) out float o_AO;
layout (location = 4) out vec4 o_SpecSH;
layout (location = 5) out vec2 o_SpecCoCg;
layout (location = 6) out float o_SpecVariance;

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
uniform sampler2D u_SpecularVarianceTexture;
uniform sampler2D u_AO;

uniform sampler2D u_SpecularSH;
uniform sampler2D u_SpecularCoCg;
uniform sampler2D u_SpecularHitDistance;
uniform sampler2DArray u_BlockPBRArray;


uniform vec2 u_Dimensions;
uniform int u_Step;
uniform int u_CurrentPass; // atrous pass iteration
uniform bool u_ShouldDetailWeight;
uniform bool DO_SPATIAL;
uniform bool u_DenoiseReflections;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform float u_ColorPhiBias = 2.0f;
uniform float u_DeltaTime;
uniform float u_Time;

const float POSITION_THRESH = 4.0f;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

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
	const float Epsilon = 0.002f;
	return v.x > 0.0f+Epsilon && v.x < 1.0f-Epsilon && v.y > 0.0f+Epsilon && v.y < 1.0f-Epsilon;
}

vec3 Saturate(vec3 x)
{
	return clamp(x, 0.0f, 1.0f);
}

float GetVarianceEstimate(out float BaseVariance)
{
#ifdef ESTIMATE_VARIANCE_BASED_ON_NEIGHBOURS
	vec2 TexelSize = 1.0f / textureSize(u_SH, 0);
	float VarianceSum = 0.0f;

	const float Kernel[3] = float[3](1.0 / 4.0, 1.0 / 8.0, 1.0 / 16.0);

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
#else 
	float x = texture(u_VarianceTexture, v_TexCoords).r;
	BaseVariance = x;
	return x;
#endif
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

int GetBlockID(vec2 txc)
{
	float id = texelFetch(u_BlockIDTexture, ivec2(txc * textureSize(u_BlockIDTexture, 0).xy), 0).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

float sqr(float x) { return x * x; }
float GetSaturation(in vec3 v) { return length(v); }

float GradientNoise()
{
	vec2 coord = gl_FragCoord.xy + mod(u_Time * 100.493850275f, 500.0f);
	float noise = fract(52.9829189f * fract(0.06711056f * coord.x + 0.00583715f * coord.y));
	return noise;
}

vec2 CalculateUV(vec3 world_pos, in vec3 normal);

void main()
{
	const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );
	const float GaussianWeights[3] = float[3]( 1.0f, 0.058141f, 0.056341f );
	const float AtrousWeights2[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);

	vec2 TexelSize = 1.0f / u_Dimensions;
	vec4 TotalColor = vec4(0.0f);

	bool FilterAO = u_Step <= 6;

	vec4 BasePosition = GetPositionAt(v_TexCoords);
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;
	vec3 BaseUtility = texture(u_Utility, v_TexCoords).xyz;
	vec4 BaseSH = texture(u_SH, v_TexCoords).xyzw;
	vec2 BaseCoCg = texture(u_CoCg, v_TexCoords).xy;
	float BaseLuminance = SHToY(BaseSH);
	float BaseVariance = 0.0f;
	float VarianceEstimate = GetVarianceEstimate(BaseVariance);
	float BaseAO = texture(u_AO, v_TexCoords).r;
	int BaseBlockID = GetBlockID(v_TexCoords);
	vec2 TexelSizeSpecSH = 1.0f / textureSize(u_SpecularSH, 0);



	// Specular Weights 
	vec4 BaseSpecularSH = texture(u_SpecularSH, v_TexCoords).xyzw;
	vec2 BaseSpecularCoCg = texture(u_SpecularCoCg, v_TexCoords).xy;
	vec2 BaseUV = CalculateUV(BasePosition.xyz, BaseNormal);
	BaseUV = clamp(BaseUV, 0.001f, 0.999f);
	float TexArrayRef = float(BlockPBRData[BaseBlockID]);
	float BaseRoughness = texture(u_BlockPBRArray, vec3(BaseUV, TexArrayRef)).r;
	float BaseSpecularLuminance = SHToY(BaseSpecularSH);







	// Start with the base inputs, one iteration of the loop can then be skipped
	vec4 TotalSH = BaseSH;
	vec2 TotalCoCg = BaseCoCg;
	float TotalWeight = 1.0f;
	float TotalVariance = BaseVariance;
	float TotalAO = BaseAO;
	float TotalAOWeight = 1.0f;
	float TotalSpecularWeight = GaussianWeights[0];//1.00000001f;
	
	float PhiColor = sqrt(max(0.0f, 1e-10 + VarianceEstimate));
	PhiColor /= max(u_ColorPhiBias, 0.4f); 

	ivec2 Jitter = ivec2((GradientNoise() - 0.5f) * float(u_Step * 1.1f));


	// Specular 

	vec4 TotalSpecSH=vec4(0.0f);
	vec2 TotalSpecCoCg=vec2(0.0f);
	TotalSpecSH += BaseSpecularSH*TotalSpecularWeight*float(BasePosition.w>0.0f);
	TotalSpecCoCg += BaseSpecularCoCg*TotalSpecularWeight*float(BasePosition.w>0.0f);

	float PassBias = mix(0.0f, 0.4f, float(u_Step<=2));
	float PassBiasLobe = mix(0.0f, 0.4f, float(u_Step<=2));

	vec2 SpecVar = vec2(0.0f);
	int SamplesTaken = 0;

	float SpecularVarianceBase = (texture(u_SpecularVarianceTexture, v_TexCoords).x);
	SpecularVarianceBase = mix(SpecularVarianceBase, 1.0f, BaseRoughness * BaseRoughness * BaseRoughness * BaseRoughness);
	float LobeFactor = mix(1.0f, 8.0f, clamp(SpecularVarianceBase / 1.1f, 0.0f, 1.0f));
	//PassBias = 0.0f;
	//PassBiasLobe = 0.0f;

	// 9 samples, with 6 passes and 1 initial pass
	for (int x = -1 ; x <= 1 ; x++)
	{
		for (int y = -2 ; y <= 2 ; y++)
		{
			vec2 SampleCoord = v_TexCoords + ((vec2(x, y) * float(u_Step)) + vec2(Jitter*0.0f)) * TexelSize;
			if (!InScreenSpace(SampleCoord)) { continue; }
			if (x == 0 && y == 0) { continue; }
			vec4 SamplePosition = GetPositionAt(SampleCoord);
			if (SamplePosition.w < 0.0f) { continue; }

			// Weights : 
			vec3 PositionDifference = abs(SamplePosition.xyz - BasePosition.xyz);
            float DistSqr = dot(PositionDifference, PositionDifference);
			vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
			float NormalWeight = pow(max(dot(BaseNormal, SampleNormal), 0.0f), 16.0f);

			

			{
				// Kernel Weights : 
				float XWeight = AtrousWeights[abs(x)];
				float YWeight = AtrousWeights[abs(y)];

				if (DistSqr < 1.4f) 
				{

					// Samples :
					vec4 SampleSH = texture(u_SH, SampleCoord).xyzw;
					vec2 SampleCoCg = texture(u_CoCg, SampleCoord).xy;
					float SampleLuma = SHToY(SampleSH);
					float SampleVariance = texture(u_VarianceTexture, SampleCoord).r;

					// :D
					float LuminosityWeight = abs(SampleLuma - BaseLuminance) / PhiColor;
					float Weight = exp(-LuminosityWeight - NormalWeight);
					Weight = max(Weight, 0.0f);

					

					Weight = (XWeight * YWeight) * Weight;
					Weight = max(Weight, 0.01f);

					TotalSH += SampleSH * Weight;
					TotalCoCg += SampleCoCg * Weight;
					TotalVariance += sqr(Weight) * SampleVariance;
					TotalWeight += Weight;

					if (FilterAO) {
						TotalAO += texture(u_AO, SampleCoord).x * Weight;
						TotalAOWeight += Weight;
					}


				}
			}



			// filter spec
			{
				// Kernel Weights : 
				float XWeight = AtrousWeights[abs(x)];
				float YWeight = AtrousWeights[abs(y)];

				const float thresh = sqrt(2.0f) * 0.9f;
				if (DistSqr < thresh && SamplePosition.w > 0.0f) 
				{

					float CurrentSpecWeight = 1.0f;
					vec2 SampleUV = CalculateUV(SamplePosition.xyz, SampleNormal);
					SampleUV = clamp(SampleUV, 0.001f, 0.999f);
					int BlockAt = GetBlockID(SampleCoord);
					
					float SampleTexArrayRef = float(BlockPBRData[BlockAt]);
					float SampleRoughness = texture(u_BlockPBRArray, vec3(SampleUV, SampleTexArrayRef)).r;
					vec4 SampleSH = texture(u_SpecularSH, SampleCoord).xyzw;
					vec2 SampleCoCg = texture(u_SpecularCoCg, SampleCoord).xy;

					float LumaAt = SHToY(SampleSH);
					float LuminanceError =  1.0f - clamp(abs(LumaAt - BaseSpecularLuminance) / 4.0f, 0.0f, 1.0f);
					float LumaWeightExponent = 1.0f;
					LumaWeightExponent = mix(0.0f, 6.0f, pow(SampleRoughness, 16.0f));
					float LuminanceWeight = pow(abs(LuminanceError), LumaWeightExponent+0.4f);

					//float LumaAt = SHToY(SampleSH);
					//float LuminanceError = 1.0f - abs(LumaAt - BaseSpecularLuminance);
					//float LumaWeightExponent = 1.0f;
					//LumaWeightExponent = mix(0.1f, 12.0f, pow(SampleRoughness, 16.0f));
					//float LuminanceWeight = pow(abs(LuminanceError), LumaWeightExponent+0.75f);
					//
					bool SampleTooRough = SampleRoughness > 0.897511f;

					LuminanceWeight = clamp(LuminanceWeight, 0.0f, 1.0f);
					NormalWeight = clamp(NormalWeight, 0.0f, 1.0f);

					// Distance weight :
					float MinRoughness = clamp(SampleRoughness, 0.0f, 1.00001f);
					float SpecularHitDistanceAt = max(texture(u_SpecularHitDistance, SampleCoord).x,0.0001f);
					vec3 ViewSpaceSample = vec3(u_View * vec4(SamplePosition.xyz, 1.0f));
					float d = length(ViewSpaceSample);
					float f = SpecularHitDistanceAt / max((SpecularHitDistanceAt + d), 1e-6f);
					float Radius = clamp(pow(mix(1.0f * MinRoughness, 1.0f, f), (1.0f - MinRoughness) * LobeFactor), 0.0f, 1.0f);
					float LobeWeight = clamp(Radius, 0.00001f, 1.0f);

					// introduces bias 
					
					CurrentSpecWeight = (XWeight * YWeight) * max(LuminanceWeight, PassBias) * NormalWeight * max(LobeWeight, PassBiasLobe);
					CurrentSpecWeight = clamp(CurrentSpecWeight, 0.0001f, 10.0f);
					TotalSpecSH += SampleSH * CurrentSpecWeight;
					TotalSpecCoCg += SampleCoCg * CurrentSpecWeight;
					TotalSpecularWeight += CurrentSpecWeight;

					float Sl = SHToY(SampleSH);
					SpecVar += vec2(Sl, Sl * Sl);
					SamplesTaken++;
				}
			}
		}
	}

	TotalSpecularWeight = max(TotalSpecularWeight, 0.001f);
	
	TotalSH /= TotalWeight;
	TotalCoCg /= TotalWeight;
	TotalVariance /= sqr(TotalWeight);
	TotalAO /= TotalAOWeight;
	TotalSpecSH /= TotalSpecularWeight;
	TotalSpecCoCg /= TotalSpecularWeight;
	
	// Output : 
	o_SH = TotalSH;
	o_CoCg = TotalCoCg;
	o_Variance = TotalVariance;
	o_AO = TotalAO;
	o_SpecSH = TotalSpecSH;
	o_SpecCoCg = TotalSpecCoCg;
	SpecVar /= max(float(SamplesTaken), 1.0f);
    float TotalSpecularVariance = max(0.0, SpecVar.y - SpecVar.x * SpecVar.x);
    o_SpecVariance = TotalSpecularVariance;

	const bool DontFilter = true;

	if (!DO_SPATIAL) { 
		o_SH = BaseSH;
		o_CoCg = BaseCoCg;
		o_Variance = BaseVariance;
		o_AO = BaseAO;
		
	}

	if (!u_DenoiseReflections) {
		o_SpecSH = BaseSpecularSH;
		o_SpecCoCg = BaseSpecularCoCg;
	}
}


bool CompareVec3(vec3 v1, vec3 v2) {
	float e = 0.0125f;
	return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
}

vec2 CalculateUV(vec3 world_pos, in vec3 normal)
{
	vec2 uv = vec2(0.0f);

	const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
	const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
	const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
	const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
	const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
	const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

    if (CompareVec3(normal, NORMAL_TOP))
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (CompareVec3(normal, NORMAL_BOTTOM))
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (CompareVec3(normal, NORMAL_RIGHT))
    {
        uv = vec2(fract(world_pos.zy));
    }

    else if (CompareVec3(normal, NORMAL_LEFT))
    {
        uv = vec2(fract(world_pos.zy));
    }
    
    else if (CompareVec3(normal, NORMAL_FRONT))
    {
        uv = vec2(fract(world_pos.xy));
    }

     else if (CompareVec3(normal, NORMAL_BACK))
    {
        uv = vec2(fract(world_pos.xy));
    }

	return uv;
}