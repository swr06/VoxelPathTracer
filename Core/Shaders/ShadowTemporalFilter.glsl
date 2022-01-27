#version 330 core

layout (location = 0) out vec4 o_Color;
layout (location = 1) out float o_Frames;


in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform sampler2D u_NormalTexture;
uniform sampler2D u_ShadowTransversals;
uniform sampler2D u_DenoisedTransversals;

uniform sampler2D u_FrameCount;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform vec3 u_PrevCameraPos;
uniform vec3 u_CurrentCameraPos;

uniform bool u_ShadowTemporal = false;
uniform bool u_ShouldFilterShadows = false;



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


bool InScreenSpace(vec2 x)
{
    return x.x < 1.0f && x.x > 0.0f && x.y < 1.0f && x.y > 0.0f;
}


bool InThresholdedScreenSpace(in vec2 v) 
{
	float b = 0.03f;
	return v.x > b && v.x < 1.0f - b && v.y > b && v.y < 1.0f - b;
}

float ManhattanDistance(vec3 p1, vec3 p2) {
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

vec3 GetShadowSpatial(float TransversalAt, float BaseDepth) {
	
	const float UnitDiagonal = sqrt(2.0f);

	if (TransversalAt < UnitDiagonal * 2.25f) {
		return vec3(1.0f);
		return texture(u_CurrentColorTexture, v_TexCoords).xyz;
	}

	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	vec3 Total = texture(u_CurrentColorTexture, v_TexCoords).xyz;
	vec3 Base = Total.xyz;
	float Weight = 1.0f;


	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;

	float Scale = 1.0f;

	Scale = TransversalAt > 8.0f ? 1.414f : 1.0f;

	for (int x = -1 ; x <= 1 ; x++) 
	{
		for (int y = -1 ; y <= 1 ; y++) 
		{
			if (x == 0 && y == 0) { continue; }

			vec2 SampleCoord = v_TexCoords + (vec2(x,y)) * TexelSize;

			if (InThresholdedScreenSpace(SampleCoord)) {
				
				float SampleDepth = texture(u_CurrentPositionTexture, SampleCoord).x;
				vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;

				if (SampleNormal == BaseNormal && abs(SampleDepth - BaseDepth) < 1.0f)
				{
					vec3 Sample = texture(u_CurrentColorTexture, SampleCoord).xyz;
					float WeightAt = clamp(1.0f - clamp(abs(Sample.x - Base.x) / 3.0f, 0.0f, 1.0f), 0.0f, 1.0f);
					WeightAt = clamp(pow(WeightAt, 6.0f), 0.0525f, 1.0f);
					Total += Sample * WeightAt;
					Weight += WeightAt;
				}
			}
		}
	}

	Total /= Weight;
	return Total;
}

vec3 clipAABB(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

vec3 ClipShadow(vec2 Reprojected) 
{
    vec3 MinColor = vec3(100.0);
	vec3 MaxColor = vec3(-100.0); 
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture,0);

	vec2 ShadowClipOffsets[5] = vec2[5](vec2(-1.0f, 0.0f), vec2(1.0f, 0.0f), vec2(0.0f, 0.0f), vec2(0.0f, -1.0f), vec2(0.0f, 1.0f));

    for(int s = 0; s < 5; s++) 
	{
        vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(ShadowClipOffsets[s]) * TexelSize).rgb; 
        MinColor = min(Sample, MinColor); 
		MaxColor = max(Sample, MaxColor); 
    }

	vec3 HistoryShadow = texture(u_PreviousColorTexture, Reprojected).xyz;
	return clipAABB(vec3(HistoryShadow.x), vec3(MinColor.x) - 0.025f, vec3(MaxColor.x) + 0.025f);
}

void main()
{
	const float Diagonal = sqrt(2.0f);

	o_Color = vec4(0.0f);

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;
		Reprojected = Reprojection(CurrentPosition.xyz);

		float ShadowTransversalAt = texture(u_ShadowTransversals, v_TexCoords).x * 100.0f;
		vec4 CurrentColor = u_ShadowTemporal ? vec4(GetShadowSpatial(ShadowTransversalAt, CurrentPosition.w), 1.0f) : texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec4 PrevColorOrig = texture(u_PreviousColorTexture, Reprojected);
		vec4 PrevColor = u_ShadowTemporal && ShadowTransversalAt < 1.414f * 3.0f ? vec4(ClipShadow(Reprojected), 1.) : PrevColorOrig;
		vec4 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyzw;

		float Bias = u_ShadowTemporal ? 0.005f : 0.01;

		if (Reprojected.x > 0.0 + Bias && Reprojected.x < 1.0 - Bias && Reprojected.y > 0.0 + Bias && Reprojected.y < 1.0 - Bias)
		{
			float d = (distance(PrevPosition.xyz, CurrentPosition.xyz));
			
			CurrentColor = clamp(CurrentColor, 0.0f, 1.0f);
			PrevColor = clamp(PrevColor, 0.0f, 1.0f);
			vec2 VelocityRejection = (v_TexCoords.xy - Reprojected.xy) * textureSize(u_CurrentColorTexture, 0).xy;

			float ClipError = abs(PrevColorOrig.x - PrevColor.x);
			float FrameIncrement = ClipError < 0.2f ? 1.0f : 0.6f;
			float FrameCountFetch = texture(u_FrameCount, Reprojected.xy).x;
			o_Frames = FrameCountFetch + FrameIncrement;

			// Linearly increasing blur factor 
			float BlendFactor = clamp((1.0f - (1.0f / o_Frames))*1.2f, 0.01f, 0.97f);
			
			BlendFactor *= clamp(exp(-length(VelocityRejection)) * 0.8f + 0.6f, 0.00000001f, 1.0f);
			
			if (d > Diagonal) {
				float DepthRejection = pow(exp(-d), 32.0f);
				BlendFactor *= clamp(DepthRejection, 0.0f, 1.0f);
			}
			
			o_Color = mix(CurrentColor, PrevColor, clamp(BlendFactor, 0.0f, 0.97f));

		}

		else 
		{
			o_Color = CurrentColor;
			o_Frames = 0.0f;
		}
	}

	else 
	{
		o_Color =  texture(u_CurrentColorTexture, v_TexCoords);
		o_Frames = 0.0f;
	}

	o_Frames = clamp(o_Frames, 0.0f, 256.0f);
}

