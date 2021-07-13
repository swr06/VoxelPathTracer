#version 330 core

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_InputTexture;
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

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
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

// Atrous weights
const float[5] Weights = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);

float ManhattanDistance(vec3 p1, vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

bool PositionDiffValid(vec3 p1, vec3 p2, float t)
{
	return abs(p1.x - p2.x) < t &&
		   abs(p1.y - p2.y) < t &&
		   abs(p1.z - p2.z) < t;
}

bool IsSky(float d)
{
	return d < 0.0f;
}

vec3 Saturate(vec3 x)
{
	return clamp(x, 0.0f, 1.0f);
}

void main()
{
	vec2 TexelSize = 1.0f / u_Dimensions;

	vec4 TotalColor = vec4(0.0f);
	float TotalWeight = 0.0f; 

	vec4 BasePosition = GetPositionAt(u_PositionTexture, v_TexCoords);
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec4 BaseColorSample = texture(u_InputTexture, v_TexCoords).rgba;
	int BaseBlock = GetBlockAt(v_TexCoords);
	float BaseSaturation = length(Saturate(BaseColorSample.xyz));
	float BaseLuminosity = GetLuminance(BaseColorSample.xyz);
	vec3 NormalizedBaseColor = normalize(Saturate(BaseColorSample.xyz));
	float BaseVariance = texture(u_VarianceTexture, v_TexCoords).r;
	bool BaseIsSky = IsSky(BasePosition.w);
	
	const float SATURATION_POW = 4.0f;
	const float DOT_POW = 1.0f;

	float VarianceWeight = 1.0f;

	//if (BaseVariance > 0.025f)
	//{
	//	VarianceWeight = 1.34f;
	//}

	for (int tx = -2 ; tx <= 2 ; tx++)
	{
		for (int ty = -2 ; ty <= 2 ; ty++)
		{
			float x = u_Step * tx;
			float y = u_Step * ty;
			float CurrentWeightX = Weights[tx + 2];
			float CurrentWeightY = Weights[ty + 2];
			vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;

			vec4 PositionAt = GetPositionAt(u_PositionTexture, SampleCoord);
			vec3 NormalAt = texture(u_NormalTexture, SampleCoord).xyz;
			int BlockAt = GetBlockAt(SampleCoord);
			vec4 SampleColor = texture(u_InputTexture, SampleCoord);
			float SamplePositionError = ManhattanDistance(PositionAt.xyz, BasePosition.xyz);

			if (SamplePositionError < 4.25f &&
				NormalAt == BaseNormal && BlockAt == BaseBlock &&
				IsSky(PositionAt.w) == BaseIsSky)
			{
				// Saturation and unit vector dot weights
				float SaturationAt = length(Saturate(SampleColor.xyz));
				vec3 NormalizedColorAt = normalize(Saturate(SampleColor.xyz));
				float SaturationWeight = pow(1.0f - abs(length(Saturate(SampleColor.xyz)) - BaseSaturation), SATURATION_POW);
				float NormalizeWeight = pow(0.5f + 0.5f * dot(NormalizedBaseColor, NormalizedColorAt), DOT_POW);
				
				// Luminosity Weights
				float CurrentLuminosity = GetLuminance(Saturate(SampleColor.xyz));
				float LuminanceError = 1.0f - abs(CurrentLuminosity - BaseLuminosity);
				float LuminanceWeight = pow(abs(LuminanceError), 8.0f);

				//float AdditionalWeight = clamp(LuminanceWeight, 0.25f, 1.5f);
				/*
				Luminosity and saturation weights are wip! 
				I still need to figure somethings out
				*/
				
				float FinalSatLumaWeight = clamp((LuminanceWeight * SaturationWeight * NormalizeWeight), 0.3f, 1.0f);
				
				// Try to reduce the overblurring by using the delta color perhaps? 
				// Delta color test : 
				vec3 DeltaColor = SampleColor.xyz - BaseColorSample.xyz;

				// Edge preservation factor :
				// This seems to work, and it works well. 
				// BUT 
				// the spatial filter's strength is reduced drastically so a lot of noise remains
				// Idea : Try making it adaptive, so you make it depend on the temporal filter as well
				// Store the temporal filter's blend factor and then have the edge preservation factor ONLY if 
				// it has blended enough frames 
				// Idea 2 : compute variance of each pixel and store in a buffer 
				// Use that as a weight (or to influence some other weight) to improve the sharpness

				float Threshold = 0.001f; 
				float x = 0.5f / (Threshold * Threshold);   
				float y = 0.3989422f / Threshold;
				float EdgePreserveFactor = exp(-dot(DeltaColor, DeltaColor) * x) * y;

				// Final Weight Sum : 
				float AdditionalWeight = 1.0f;
				if (BaseVariance < 0.005f && u_ShouldDetailWeight) {
					AdditionalWeight = clamp(EdgePreserveFactor, 0.08500f, 1.0f);
				}

				float DistanceWeight = abs(4.25f - SamplePositionError);
				DistanceWeight = pow(DistanceWeight, 5.0f); // Strong attenuation

				float Weight = CurrentWeightX * CurrentWeightY *
							   DistanceWeight *
							   clamp(AdditionalWeight, 0.0f, 1.0f) *
							   VarianceWeight;

				TotalColor += SampleColor * Weight;
				TotalWeight += Weight;
			}
		}
	}

	o_SpatialResult = TotalColor / max(TotalWeight, 0.01f);
	//o_SpatialResult = texture(u_InputTexture, v_TexCoords).rgba;
}