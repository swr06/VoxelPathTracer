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
uniform sampler2D u_GBufferNormals;
uniform sampler2D u_GBufferPBR;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform vec2 u_Dimensions;
uniform int u_Step;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

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


// Basic clamp firefly reject
vec4 FireflyReject(vec4 Col)
{
	return vec4(clamp(Col.xyz, 0.0f, 1.0f + 0.8f), Col.w);
}

int GetBlockID(vec2 txc)
{
	float id = texelFetch(u_BlockIDTex, ivec2(txc * textureSize(u_BlockIDTex, 0).xy), 0).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

//bool SampleNormalMappedAt(vec3 WorldPos, out vec3 N) {
//	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
//	ProjectedPosition.xyz /= ProjectedPosition.w;
//	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;
//	float b = 0.02f;
//	bool v = ProjectedPosition.x > b && ProjectedPosition.x < 1.0f - b && ProjectedPosition.y > b && ProjectedPosition.y < 1.0f - b;
//	N = texture(u_NormalMappedTexture, ProjectedPosition.xy).xyz;
//	return v;
//}
//


void main()
{
	vec4 BlurredSH = vec4(0.0f);
	vec2 BlurredCoCg = vec2(0.0f);

	vec3 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;
	int BaseBlockID = GetBlockID(v_TexCoords);

	vec4 BaseSH = texture(u_InputTexture, v_TexCoords).xyzw;
	float BaseLuminance = SHToY(BaseSH);


	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	float TexArrayRef = float(BlockPBRData[BaseBlockID]);
	float RoughnessAt = texture(u_GBufferPBR, v_TexCoords).r;

	for (int s = 0 ; s < GAUSS_KERNEL; s++)
	{
		int Sample = GaussianOffsets[s]; // todo : use u_Step here!
		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (Sample * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (Sample * TexelSize));
		
		// Solves clamp issues : 
		float bias = 0.01f;
		if (SampleCoord.x > 0.0f + bias && SampleCoord.x < 1.0f - bias && SampleCoord.y > 0.0f + bias && SampleCoord.y < 1.0f - bias) 
		{

			vec3 SamplePosition = GetPositionAt(u_PositionTexture, SampleCoord).xyz;
			vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;

			vec3 PositionDifference = abs(SamplePosition - BasePosition);
			float PositionError = dot(PositionDifference, PositionDifference);

			int BlockAt = GetBlockID(SampleCoord);

			if (PositionError > 0.75f || SampleNormal != BaseNormal || BlockAt != BaseBlockID) { 
				continue;
			}

			vec4 SampleSH = texture(u_InputTexture, SampleCoord).xyzw;
			vec2 SampleCoCg = texture(u_InputCoCgTexture, SampleCoord).rg;

			// Luminosity weights
			float LumaAt = SHToY(SampleSH);
			float LuminanceError = 1.0f - abs(LumaAt - BaseLuminance);


			//float LumaTolerance = mix(8.0f, 0.0f, clamp(RoughnessAt * 2.0f, 0.0f, 1.0f));
			//LumaTolerance = clamp(LumaTolerance, 0.25f, 8.0f);


			float LumaWeightExponent = 1.0f;
			LumaWeightExponent = mix(0.1f, 8.0f, RoughnessAt*RoughnessAt);
			float LuminanceWeight = pow(abs(LuminanceError), LumaWeightExponent);
			float CurrentKernelWeight = GaussianWeightsNormalized[s];

			float CurrentWeight = 1.0f;


			CurrentWeight = CurrentKernelWeight * CurrentWeight;
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

		// Preserve a few more details in smoother materials
		if (RoughnessAt <= 0.125f)
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

