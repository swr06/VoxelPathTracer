// Shadow/Diffuse/Old Reflection Temporal filter ->

#version 330 core

layout (location = 0) out vec4 o_Color;
layout (location = 1) out vec2 o_CoCg; // Optional output 

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

// (Actually the CoCg part, the chrominance is stored in the first band of the spherical harmonic[SHy]) 
// idk wtf happened to me when naming this 
uniform sampler2D u_PreviousSH;
uniform sampler2D u_CurrentSH;

uniform sampler2D u_NormalTexture;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_MinimumMix = 0.25f;
uniform float u_MaximumMix = 0.975f;
uniform int u_TemporalQuality = 1; // 0, 1, 2

uniform vec3 u_PrevCameraPos;
uniform vec3 u_CurrentCameraPos;

uniform bool u_ReflectionTemporal = false;
uniform bool u_DiffuseTemporal = false;
uniform bool u_ShadowTemporal = false;
uniform bool u_ShouldFilterShadows = false;
uniform float u_ClampBias = 0.025f;


vec2 Dimensions;

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

vec4 ClipColor(vec4 aabbMin, vec4 aabbMax, vec4 prevColor) 
{
    vec4 pClip = (aabbMax + aabbMin) / 2;
    vec4 eClip = (aabbMax - aabbMin) / 2;
    vec4 vClip = prevColor - pClip;
    vec4 vUnit = vClip / eClip;
    vec4 aUnit = abs(vUnit);
    float divisor = max(aUnit.x, max(aUnit.y, aUnit.z));

    if (divisor > 1)
	{
        return pClip + vClip / divisor;
    }

    return prevColor;
}

vec4 GetClampedColor(vec2 reprojected, in vec3 worldpos)
{
	int quality = clamp(u_TemporalQuality, 0, 2);
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	ivec2 Coord = ivec2(v_TexCoords * Dimensions);
	
	if (quality == 1)
	{

		vec4 minclr = vec4(10000.0f); 
		vec4 maxclr = vec4(-10000.0f); 

		const int SampleSize = 1;

		for(int x = -SampleSize; x <= SampleSize; x++) 
		{
			for(int y = -SampleSize; y <= SampleSize; y++) 
			{
				vec4 Fetch = texelFetch(u_CurrentColorTexture, Coord + ivec2(x,y), 0); 
				minclr = min(minclr, Fetch); 
				maxclr = max(maxclr, Fetch); 
			}
		}

		minclr -= u_ClampBias; 
		maxclr += u_ClampBias; 
		
		return clamp(texture(u_PreviousColorTexture, reprojected), minclr, maxclr); 
	}

	else if (quality == 2)
	{
		vec2 TexelSize = 1.0f / textureSize(u_PreviousColorTexture, 0).xy;
		
		vec2 BestOffset = vec2(0.0f, 0.0f);
		float BestDiff = 10000.0f;

		vec4 minclr = vec4(10000.0f);
		vec4 maxclr = vec4(-10000.0f);

		const int BoxSampleSize = 2;

		for(int x = -BoxSampleSize; x <= BoxSampleSize; x++) 
		{
			for(int y = -BoxSampleSize; y <= BoxSampleSize; y++) 
			{
				vec4 SampledPosition = GetPositionAt(u_PreviousFramePositionTexture, reprojected + (vec2(x, y) * TexelSize)).rgba;
				vec4 Fetch = texelFetch(u_CurrentColorTexture, Coord + ivec2(x,y), 0); 

				minclr = min(minclr, Fetch.xyzw); 
				maxclr = max(maxclr, Fetch.xyzw); 

				if (SampledPosition.w > 0.0f)
				{
					float Diff = abs(distance(worldpos, SampledPosition.xyz));

					if (Diff < BestDiff)
					{
						BestDiff = Diff;
						BestOffset = vec2(x, y);
					}
				}	
			}
		}

		minclr -= u_ClampBias; 
		maxclr += u_ClampBias; 

		vec4 FinalColor = texture(u_PreviousColorTexture, reprojected + (BestOffset * TexelSize)).xyzw;
		return clamp(FinalColor, minclr, maxclr);
	}

	else 
	{
		return texture(u_PreviousColorTexture, reprojected).rgba;
	}
}

bool InThresholdedScreenSpace(in vec2 v) 
{
	float b = 0.03f;
	return v.x > b && v.x < 1.0f - b && v.y > b && v.y < 1.0f - b;
}

float ManhattanDistance(vec3 p1, vec3 p2) {
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

vec4 GetShadowSpatial() {
	
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	vec4 Total = texture(u_CurrentColorTexture, v_TexCoords);
	float Samples = 1.0f;

	vec3 BasePosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;

	float Scale = 1.5f;

	for (int x = -1 ; x <= 1 ; x++) 
	{
		for (int y = -2 ; y <= 2 ; y++) 
		{
			if (x == 0 && y == 0) { continue; }

			vec2 SampleCoord = v_TexCoords + (vec2(x,y)*Scale) * TexelSize;

			if (InThresholdedScreenSpace(SampleCoord)) {
				
				vec3 SamplePosition = GetPositionAt(u_CurrentPositionTexture, SampleCoord).xyz;
				vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;

				if (SampleNormal == BaseNormal && ManhattanDistance(SamplePosition, BasePosition) < 1.414141414f)
				{
					Total += texture(u_CurrentColorTexture, v_TexCoords).xyzw;
					Samples += 1.0f;
				}
			}
		}
	}

	Total /= Samples;
	return Total;
}


// Variance clip -> aabb 
// Calculate standard deviation with average color and clip to an aabb
vec3 clipToAABB(in vec3 cOld, in vec3 cNew, in vec3 centre, in vec3 halfSize)
{
    if (all(lessThanEqual(abs(cOld - centre), halfSize))) {
        return cOld;
    }
    
    vec3 dir = (cNew - cOld);
    vec3 near = centre - sign(dir) * halfSize;
    vec3 tAll = (near - cOld) / dir;
    float t = 1e20;
    for (int i = 0; i < 3; i++) {
        if (tAll[i] >= 0.0 && tAll[i] < t) {
            t = tAll[i];
        }
    }
    
    if (t >= 1e20) {
		return cOld;
    }
    return cOld + dir * t;
}


// rgb -> ycocg and 
// ycocg -> rgb conversions 
vec3 rgb2ycocg(in vec3 rgb)
{
    float co = rgb.r - rgb.b;
    float t = rgb.b + co / 2.0;
    float cg = rgb.g - t;
    float y = t + cg / 2.0;
    return vec3(y, co, cg);
}


vec3 ycocg2rgb(in vec3 ycocg)
{
    float t = ycocg.r - ycocg.b / 2.0;
    float g = ycocg.b + t;
    float b = t - ycocg.g / 2.0;
    float r = ycocg.g + b;
    return vec3(r, g, b);
}
////

// AABB color clipping 
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


// reduces ghosting 
// spatial filter does carry quite a bit 
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

	 // Bias to avoid noticable clipping, while at the same time reducing ghosting
	return clipAABB(HistoryShadow, MinColor-0.05f, MaxColor+0.05f);
}

// Clips based on variance 
vec3 VarianceClip(vec2 Reproj, vec3 CenterCurrent) 
{
	vec2 Offsets[4] = vec2[4](vec2(-1.0,  0.0), 
								vec2( 1.0,  0.0), 
								vec2( 0.0, -1.0), 
								vec2( 0.0,  1.0));
	
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	vec3 C = texture(u_PreviousColorTexture, Reproj).xyz;
	vec3 Averaged = rgb2ycocg(CenterCurrent.rgb);
    vec3 StandardDeviation = Averaged * Averaged;

    for (int i = 0; i < 4; i++) 
	{
        vec3 Sample = rgb2ycocg(texture(u_CurrentColorTexture, (Reproj + Offsets[i] * TexelSize)).rgb);
        Averaged += Sample;
        StandardDeviation += Sample * Sample;
    }
	
    Averaged /= 5.0f;
    StandardDeviation = sqrt(StandardDeviation / 5.0f - Averaged * Averaged);
	vec3 ClippedColor = C;
	ClippedColor = ycocg2rgb(clipToAABB(rgb2ycocg(ClippedColor), rgb2ycocg(CenterCurrent.rgb), Averaged, StandardDeviation));
	return ClippedColor;
}


// I know i use a shitton of uniform branches but they are nearly free so i dont really care 
void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;
		Reprojected = Reprojection(CurrentPosition.xyz);

		vec4 CurrentColor = u_ShadowTemporal ? GetShadowSpatial() : texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec4 PrevColor = u_ShadowTemporal ? vec4(ClipShadow(Reprojected),1.0f) : texture(u_PreviousColorTexture, Reprojected);
		vec3 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyz;

		float Bias = u_ShadowTemporal ? 0.006f : 0.01;

		if (Reprojected.x > 0.0 + Bias && Reprojected.x < 1.0 - Bias && Reprojected.y > 0.0 + Bias && Reprojected.y < 1.0 - Bias)
		{
			float d = abs(distance(PrevPosition, CurrentPosition.xyz));
			float t = u_ShadowTemporal ? 0.5f - 0.01f : 1.1f;

			if (d > t) 
			{
				o_Color = CurrentColor;

				if (u_DiffuseTemporal) {
					vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
					o_CoCg = CurrentSH;
				}

				return;
			}

			float BlendFactor = d;
			BlendFactor = exp(-BlendFactor);
			BlendFactor = clamp(BlendFactor, clamp(u_MinimumMix, 0.01f, 0.9f), clamp(u_MaximumMix, 0.1f, 0.98f));
			
			if (u_ShadowTemporal) {  
				CurrentColor = clamp(CurrentColor + 0.005f, 0.0f, 1.0f);  // Bias
				PrevColor = clamp(PrevColor + 0.005f, 0.0f, 1.0f);  // Bias

				if (!u_ShouldFilterShadows) {
					BlendFactor = 0.0f;
				}
			}

			o_Color = mix(CurrentColor, PrevColor, BlendFactor);

			if (u_ShadowTemporal) {  
				o_Color = clamp(o_Color, 0.0f, 1.0f);  // Bias
			}

			if (u_DiffuseTemporal) {
				vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
				vec2 PrevSH = texture(u_PreviousSH, Reprojected.xy).xy;
				o_CoCg = mix(CurrentSH, PrevSH, BlendFactor);
			}
		}

		else 
		{
			o_Color = CurrentColor;

			if (u_DiffuseTemporal) {
				vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
				o_CoCg = CurrentSH;
			}
		}
	}

	else 
	{
		o_Color = u_ShadowTemporal ? GetShadowSpatial() : texture(u_CurrentColorTexture, v_TexCoords);
		if (u_DiffuseTemporal) {
				vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
				o_CoCg = CurrentSH;
		}
	}
}

