// Gaussian reflection denoiser
// This denoiser has specialized weights so that it doesn't overblur anything 

// "Welcome to magic number land, we hope you enjoy your stay."


#version 430 core

//#define NORMAL_MAP_KERNEL_WEIGHT

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_BlockIDTex;
uniform sampler2D u_SpecularHitData;

uniform sampler2D u_GBufferNormals;
uniform sampler2D u_GBufferPBR;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform bool u_RoughnessBias; 
uniform bool u_NormalMapAware; 
uniform bool u_HandleLobeDeviation; 
uniform bool u_DeriveFromDiffuseSH; 
uniform vec2 u_Dimensions;
uniform int u_Step;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform float u_Time;

uniform int u_ReflectionDenoisingRadiusBias;

uniform float u_NormalMapWeightStrength;
uniform float u_ReflectionDenoiserScale;

uniform float u_ResolutionScale;

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
	// Calculate standard deviation 
	// Mean is the mean of a box kernel on a 3x3 or 5x5 area
	// While calculating the mean, average the luminance squared as well to calculate variance.

    vec3 StandardDeviation = vec3(sqrt(max(0.00001f, VarianceEstimate))); 
    vec3 Threshold = 0.1f + Mean + StandardDeviation * 8.0;
    vec3 ErrorEstimate = vec3(max(vec3(0.0f), Radiance - Threshold));

	// dim based on variance error
    return clamp(Radiance - ErrorEstimate, 0.0f, 16.0f); 
}

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

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}


void main()
{
	const float Diagonal = sqrt(2.0f);
	
	vec4 FilteredColor = vec4(0.0f);
	vec2 BlurredCoCg = vec2(0.0f);

	vec4 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords).xyzw;
	
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;

	int BaseBlockID = GetBlockID(v_TexCoords);

	bool BaseIsSky = BasePosition.w < 0.0f;
	vec4 BaseColor = texture(u_InputTexture, v_TexCoords).xyzw;
	float BaseLuminance = GetLuminance(BaseColor.xyz);

	float TotalWeight = 0.0f;
	float TexelSize = u_Dir ? 1.0f / u_Dimensions.x : 1.0f / u_Dimensions.y;
	float TexArrayRef = float(BlockPBRData[BaseBlockID]);

	vec3 SampledPBR = texture(u_GBufferPBR, v_TexCoords).xyz;
	float RoughnessAt = SampledPBR.x;
	RoughnessAt += 0.001f;
	float RawRoughness = RoughnessAt;
	RoughnessAt *= mix(1.0f, 0.91f, float(u_RoughnessBias));

	vec3 NormalMappedBase = texture(u_GBufferNormals, v_TexCoords).xyz;

	float HitDistanceFetch = texture(u_SpecularHitData, v_TexCoords).x + 0.0001f;
	float RawHitDistance = HitDistanceFetch;

	const float LobeDistanceCurveBias = 1.0f / 1.3f;


	// Tweak transversal ->
	if (HitDistanceFetch < 0.001f) {
		HitDistanceFetch = 1.75f;
	}

	else {
		HitDistanceFetch = pow(HitDistanceFetch, LobeDistanceCurveBias);
	}

	// Handle roughness lobe deviation

	if (u_HandleLobeDeviation) {
		if (RoughnessAt <= 0.25f + 0.05f) {
			HitDistanceFetch = clamp(HitDistanceFetch, 0.0f, 8.0f);
		}

		if (RoughnessAt <= 0.2) {
			HitDistanceFetch = clamp(HitDistanceFetch, 0.0f, 5.5f);
		}
	}

	// Distance weight ->
	float SpecularHitDistance = max(HitDistanceFetch,0.02f) / 1.3f;

	// If the normal map weight is enabled, we can get decent results by tweaking the view space position 
	NormalMappedBase.z *= 1.5f;
	NormalMappedBase.x *= 4.5f;


	vec3 ViewSpaceBase = vec3(u_View * vec4(BasePosition.xyz + (NormalMappedBase * 0.75f), 1.0f));
	float ViewLength = length(ViewSpaceBase);
	float ViewLengthWeight = 0.001f + ViewLength;


	// View length weight ->

	if (RoughnessAt > 0.135f) {
		ViewLengthWeight = max(ViewLengthWeight, 0.750);
	}

	else {
		ViewLengthWeight = max(ViewLengthWeight, 3.0f);
	}

	if (RoughnessAt < 0.125f) {
		ViewLengthWeight = clamp(ViewLengthWeight, 0.000001f, 6.0f);
	}

	else if (RoughnessAt < 0.25f) {
		ViewLengthWeight = clamp(ViewLengthWeight, 0.000001f, 8.0f+1.0f);
	}

	else if (RoughnessAt < 0.5f) {
		ViewLengthWeight = clamp(ViewLengthWeight, 0.000001f, 16.0f);
	}

	else if (RoughnessAt < 0.75f) {
		ViewLengthWeight = clamp(ViewLengthWeight, 0.000001f, 24.0f);
	}

	else {
		ViewLengthWeight = clamp(ViewLengthWeight, 0.000001f, 32.0f);
	}

	float TransversalContrib = SpecularHitDistance / max((SpecularHitDistance + ViewLengthWeight), 0.00001f);
	float TransversalContrib_ = TransversalContrib;

	if (RawHitDistance < 3.5f && RawHitDistance > 0.0000001f && RoughnessAt <= 0.625f) {
		TransversalContrib = pow(TransversalContrib, 2.5f);
	}

	if (RawHitDistance < 4.0f && RawHitDistance > 0.0000001f && RoughnessAt <= 0.525f) {
		TransversalContrib = pow(TransversalContrib, 1.35f);
	}

	if (RawHitDistance < 4.0f && RawHitDistance > 0.0000001f && RoughnessAt <= 0.35f) {
		TransversalContrib = pow(TransversalContrib, 1.4f);
	}


	float Radius = clamp(pow(mix(1.0f * RoughnessAt, 1.0f, TransversalContrib), pow((1.0f-RoughnessAt),1.0/1.4f)*5.0f), 0.0f, 1.0f);
	
	// Calculate gaussian radius ->
	int EffectiveRadius = int(floor(Radius * 15.0f));
	EffectiveRadius = clamp(EffectiveRadius, 1, 15);

	// reduce noise on too rough objects 
	bool BaseTooRough = RoughnessAt > 0.897511f;
	EffectiveRadius = BaseTooRough ? 15 : EffectiveRadius;

	// Basic jitter to reduce 2 pass gaussian artifacts 
	int Jitter = int((GradientNoise() - 0.5f) * 1.25f);
	
	// Sample scale based on resolution scale
	float Scale = 1.0f;
	Scale = mix(1.0f, 2.0f, clamp(u_ResolutionScale, 0.0000001f, 1.0f)) + 0.5f;
	
	int RadiusBias = 0;
	
	if (RoughnessAt > 0.4f) {
		RadiusBias = 1;
	}

	if (RoughnessAt > 0.5f) {
		RadiusBias += 1;
	}

	if (SampledPBR.y > 0.1f && RoughnessAt > 0.45f) {
		RadiusBias += 1;
	}
	
	EffectiveRadius = clamp(EffectiveRadius + RadiusBias + u_ReflectionDenoisingRadiusBias,1,15);
	
	
	


	//Scale = mix(1.0f, 2.0f, (clamp(pow(mix(1.0f * RoughnessAt, 1.0f, TransversalContrib_), pow((1.0f-RoughnessAt),1.0/2.4f)*8.0f), 0.0f, 1.0f) * 1.75f)+0.4f);

	// Since we derived this pixel from the diffuse spherical harmonic (which is already denoised)
	// We can improve performance by limiting the spatial filter's radius and using a higher scale
	if (u_DeriveFromDiffuseSH && RawRoughness >= 0.865f) {
		EffectiveRadius = 4; // 8 spatial samples 
		Scale *= 1.25f;
	}

	Scale *= u_ReflectionDenoiserScale;

	EffectiveRadius = clamp(EffectiveRadius,1,15);

	// -- Gaussian filter -- 

	for (int Sample = -EffectiveRadius ; Sample <= EffectiveRadius; Sample++)
	{
		float SampleOffset = Sample;
		vec2 SampleCoord = u_Dir ? vec2(v_TexCoords.x + (SampleOffset * Scale * TexelSize), v_TexCoords.y) : vec2(v_TexCoords.x, v_TexCoords.y + (SampleOffset * Scale * TexelSize));
		vec2 Mask = u_Dir ? vec2(1.0f, 0.0f) : vec2(0.0, 1.0f);
		
		float bias = 0.01f;
		if (SampleCoord.x > 0.0f + bias && SampleCoord.x < 1.0f - bias && SampleCoord.y > 0.0f + bias && SampleCoord.y < 1.0f - bias) 
		{

			//vec4 SamplePosition = GetPositionAt(u_PositionTexture, SampleCoord).xyzw;
			float SampleDepth = texture(u_PositionTexture, SampleCoord).x;
			int BlockAt = GetBlockID(SampleCoord);
			bool BlockValidity = BlockAt == BaseBlockID;
			bool SampleIsSky = SampleDepth < 0.0f;
			
			if (SampleIsSky != BaseIsSky) {
				continue;
			}

			bool KillDepthWeight = false;

			if (BasePosition.w < 42.0f) {
			
				float PositionTolerance = BlockValidity ? 0.925f : 1.414f;

				float PositionError = abs(SampleDepth-BasePosition.w);
			
				
				if (PositionError > PositionTolerance || SampleDepth < 0.0f) { 
					continue;
				}

				KillDepthWeight = true;
			}
			
			

			float DepthDifference = abs(SampleDepth-BasePosition.w);
			float DepthWeight = mix(clamp(pow(exp(-max(DepthDifference, 0.00001f)), 3.0f), 0.0001f, 1.0f),1.0f,float(KillDepthWeight));
			
			vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;

			if (SampleNormal != BaseNormal) {
				continue;
			}

			//float NormalWeight = pow(abs(dot(SampleNormal, BaseNormal)), 12.0f);
			float SampleRoughness = texture(u_GBufferPBR, SampleCoord).r;

			

			vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;

			// Luminosity weights
			//float LumaAt = SHToY(SampleSH);
			//float LuminanceError = 1.0f - abs(LumaAt - BaseLuminance);
			//float LumaWeightExponent = 1.0f;
			//LumaWeightExponent = mix(0.1f, 8.0f, pow(SampleRoughness, 16.0f));
			//float LuminanceWeight = pow(abs(LuminanceError), LumaWeightExponent+0.75f);

			float LumaAt = GetLuminance(SampleColor.xyz);
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

			CurrentWeight *= DepthWeight;



			if (u_NormalMapAware && !SampleTooRough) {
				vec3 NormalMapAt = texture(u_GBufferNormals, SampleCoord).xyz;
				CurrentWeight *= pow(clamp(dot(NormalMapAt,NormalMappedBase), 0.00001f, 1.0f), 64.0f * u_NormalMapWeightStrength);
			}


			//CurrentWeight *= clamp(NormalWeight, 0.0f, 1.0f);
			CurrentWeight = clamp(CurrentWeight,0.01,1.0f);
			CurrentWeight = CurrentKernelWeight * CurrentWeight;
			CurrentWeight = clamp(CurrentWeight,0.01,1.0f);
			FilteredColor += SampleColor * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	vec4 BaseSampled = texture(u_InputTexture, v_TexCoords).rgba;

	bool DoSpatial = true;

	if (TotalWeight > 0.001f && DoSpatial) { 
		FilteredColor = FilteredColor / TotalWeight;
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

		o_SpatialResult = mix(BaseSampled, FilteredColor, m);
	}

	else {
		o_SpatialResult = BaseSampled;
	}
}
