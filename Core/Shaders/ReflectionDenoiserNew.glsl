// Gaussian reflection denoiser
// This denoiser has specialized weights so that it doesn't overblur anything 

#version 430 core

layout (location = 0) out vec4 o_SpatialResult;
layout (location = 1) out vec2 o_CoCg;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_InputTexture;
uniform sampler2D u_InputCoCgTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_BlockIDTex;
uniform sampler2D u_NormalMappedTexture;
uniform sampler2D u_SpecularHitData;
uniform sampler2DArray u_BlockPBRTexArray;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform vec2 u_Dimensions;
uniform int u_Step;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform float u_Time;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

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
	0.058754, // 16 idx
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

float SHToY(vec4 shY)
{
    return max(0, 3.544905f * shY.w);
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

vec3 VarianceFireflyRejection(vec3 Radiance, float VarianceEstimate, vec3 Mean)
{
    vec3 StandardDeviation = vec3(sqrt(max(0.00001f, VarianceEstimate))); // Calculate standard deviation 
    vec3 Threshold = 0.1f + Mean + StandardDeviation * 8.0;
    vec3 ErrorEstimate = vec3(max(vec3(0.0f), Radiance - Threshold));
    return clamp(Radiance - ErrorEstimate, 0.0f, 16.0f); // -> Dim based on variance error 
}

vec2 CalculateUV(vec3 world_pos, in vec3 normal);

int GetBlockID(vec2 txc)
{
	float id = texelFetch(u_BlockIDTex, ivec2(txc * textureSize(u_BlockIDTex, 0).xy), 0).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

float GradientNoise()
{
	vec2 coord = gl_FragCoord.xy + mod(u_Time * 100.493850275f, 500.0f);
	float noise = fract(52.9829189f * fract(0.06711056f * coord.x + 0.00583715f * coord.y));
	return noise;
}

void main()
{
	vec4 BlurredSH = vec4(0.0f);
	vec2 BlurredCoCg = vec2(0.0f);

	vec4 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords).xyzw;
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;
	int BaseBlockID = GetBlockID(v_TexCoords);
	bool BaseIsSky = BasePosition.w < 0.0f;
	vec4 BaseSH = texture(u_InputTexture, v_TexCoords).xyzw;
	float BaseLuminance = SHToY(BaseSH);

	vec2 BaseUV = CalculateUV(BasePosition.xyz, BaseNormal);
	BaseUV = clamp(BaseUV, 0.001f, 0.999f);
	BaseUV.y = 1.0f - BaseUV.y;

	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	float TexArrayRef = float(BlockPBRData[BaseBlockID]);
	float RoughnessAt = texture(u_BlockPBRTexArray, vec3(BaseUV, TexArrayRef)).r;
	float SpecularHitDistance = max(texture(u_SpecularHitData, v_TexCoords).x,0.0f)/1.3f;
	vec3 ViewSpaceBase = vec3(u_View * vec4(BasePosition.xyz, 1.0f));
	float d = length(ViewSpaceBase);
	float f = SpecularHitDistance / max((SpecularHitDistance + d), 1e-6f);
	float Radius = clamp(pow(mix(1.0f * RoughnessAt, 1.0f, f), (1.0f-RoughnessAt)*4.0f), 0.0f, 1.0f);
	int EffectiveRadius = int(floor(Radius * 12.0f));
	EffectiveRadius = clamp(EffectiveRadius, 1, 12+2);

	// reduce noise on too rough objects 
	bool BaseTooRough = RoughnessAt > 0.897511f;
	EffectiveRadius = BaseTooRough ? 12 : EffectiveRadius;
	int Jitter = int((GradientNoise() - 0.5f) * 1.25f);
	float Scale = 1.225f;

	EffectiveRadius = clamp(EffectiveRadius,1,14);

	for (int Sample = -EffectiveRadius ; Sample <= EffectiveRadius; Sample++)
	{
		float SampleOffset = Jitter + Sample;
		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (SampleOffset * Scale * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (SampleOffset * TexelSize));
		
		float bias = 0.01f;
		if (SampleCoord.x > 0.0f + bias && SampleCoord.x < 1.0f - bias && SampleCoord.y > 0.0f + bias && SampleCoord.y < 1.0f - bias) 
		{

			vec4 SamplePosition = GetPositionAt(u_PositionTexture, SampleCoord).xyzw;
			vec3 PositionDifference = abs(SamplePosition.xyz - BasePosition.xyz);
			float PositionError = dot(PositionDifference, PositionDifference);
			if (PositionError > 1.069420f || SamplePosition.w < 0.0f) { 
				continue;
			}
			
			vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
			float NormalWeight = pow(abs(dot(SampleNormal, BaseNormal)), 12.0f);
			int BlockAt = GetBlockID(SampleCoord);
			vec2 SampleUV = CalculateUV(SamplePosition.xyz, SampleNormal);
			SampleUV = clamp(SampleUV, 0.001f, 0.999f);
			float SampleTexArrayRef = float(BlockPBRData[BlockAt]);
			float SampleRoughness = texture(u_BlockPBRTexArray, vec3(SampleUV, SampleTexArrayRef)).r;

			vec4 SampleSH = texture(u_InputTexture, SampleCoord).xyzw;
			vec2 SampleCoCg = texture(u_InputCoCgTexture, SampleCoord).rg;

			// Luminosity weights
			//float LumaAt = SHToY(SampleSH);
			//float LuminanceError = 1.0f - abs(LumaAt - BaseLuminance);
			//float LumaWeightExponent = 1.0f;
			//LumaWeightExponent = mix(0.1f, 8.0f, pow(SampleRoughness, 16.0f));
			//float LuminanceWeight = pow(abs(LuminanceError), LumaWeightExponent+0.75f);

			float LumaAt = SHToY(SampleSH);
			float LuminanceError =  1.0f - clamp(abs(LumaAt - BaseLuminance) / 4.0f, 0.0f, 1.0f);
			float LumaWeightExponent = 1.0f;
			LumaWeightExponent = mix(0.01f, 8.0f, pow(SampleRoughness, 16.0f));
			float LuminanceWeight = pow(abs(LuminanceError), LumaWeightExponent + 0.5525250f);
			float CurrentKernelWeight = GaussianWeightsNormalized[clamp(16+Sample,0,32)];

			float CurrentWeight = 1.0f;

			// dont weight by luminance error if the material smoothness < threshold 
			bool SampleTooRough = SampleRoughness > 0.897511f;
			if (!SampleTooRough) {
				CurrentWeight *= clamp(LuminanceWeight, 0.0f, 1.0f);
			}

			CurrentWeight *= clamp(NormalWeight, 0.0f, 1.0f);
			CurrentWeight = clamp(CurrentWeight,0.01,1.0f);
			CurrentWeight = CurrentKernelWeight * CurrentWeight;
			CurrentWeight = clamp(CurrentWeight,0.01,1.0f);
			BlurredSH += SampleSH * CurrentWeight;
			BlurredCoCg += SampleCoCg * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	vec4 BaseSampled = texture(u_InputTexture, v_TexCoords).rgba;
	vec2 BaseSampledCoCg = texture(u_InputCoCgTexture, v_TexCoords).rg;

	bool DoSpatial = true;

	if (TotalWeight > 0.001f && DoSpatial) { 
		BlurredSH = BlurredSH / TotalWeight;
		BlurredCoCg = BlurredCoCg / TotalWeight;

		float m = 1.0f;
		
		// Preserve a few more details in "too" smooth materials
		if (RoughnessAt <= 0.10550f)
		{
			m = RoughnessAt * 16.0f;
			m = 1.0f - m;
			m = m * m * m * m;
			m = clamp(m, 0.01f, 0.999f); 
		}  

		o_SpatialResult = mix(BaseSampled, BlurredSH, m);
		o_CoCg = mix(BaseSampledCoCg, BlurredCoCg, m);
	}

	else {
		o_SpatialResult = BaseSampled;
		o_CoCg = BaseSampledCoCg;
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