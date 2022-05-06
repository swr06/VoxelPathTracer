// Todo : Handle sky motion vectors and specular motion vectors

#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

// Current frame textures
uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_PositionTexture;

// History
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousPositionTexture;

// Matrices 
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InversePrevProjection;
uniform mat4 u_InversePrevView;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_View;
uniform mat4 u_Projection;

// TAA flag 
uniform bool u_Enabled;

// true if the world was modified this frame 
uniform bool u_BlockModified;

// History positions 
uniform vec3 u_CameraHistory[2];

// Depth rejection flags/modifiers 
uniform bool u_DepthWeight;
uniform float u_DepthExponentMultiplier;


// Samplers 
vec2 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;
	return ProjectedPosition.xy;
}

vec3 GetRayDirection(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SampleBasePosition(sampler2D pos_tex)
{
	float Dist = 1./texture(pos_tex, v_TexCoords).r;
	return vec4(v_RayOrigin + normalize(v_RayDirection) * Dist, Dist);
}

vec3 GetPREVRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InversePrevProjection * clip), -1.0, 0.0);
	return vec3(u_InversePrevView * eye);
}

vec4 GetPREVPositionAt(vec2 txc)
{
	float Dist = 1./texture(u_PreviousPositionTexture, txc).r;
	return vec4(u_InversePrevView[3].xyz + normalize(GetPREVRayDirectionAt(txc)) * Dist, Dist);
}

// -- 

// AABB Clipping
vec3 AABBClip(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

// Percieved luma
float Luminance(vec3 RGB )
{
    return dot(RGB, vec3(0.2126f, 0.7152f, 0.0722f));
}

// Reinhard transform
const float ReinhardExp = 1.000000f;
vec3 Reinhard(vec3 RGB )
{
	RGB *= ReinhardExp;
    return vec3(RGB) / (vec3(1.0f) + Luminance(RGB));
}

// Inverse reinhard transform
vec3 InverseReinhard(vec3 RGB)
{
    return (RGB / (vec3(1.0f) - Luminance(RGB))) / ReinhardExp;
}

// Gets the nearest depth sample in a 3x3 area 
vec2 GetNearestDepths3x4(vec2 Reprojected, out bool SkyEdge) 
{
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0).xy;

	float MinCurrentDepth = 1000.0f;
	float MinEffectiveCurrentDepth = 2000.0f;
	float MinPreviousDepth = 1000.0f;
	float MinEffectivePreviousDepth = 2000.0f;
	
	// Sky edge flag 
	SkyEdge = false;

	// The gbuffer depth is stored differently, hence, it needs slightly more operations to get the min/max depths.
	for (int x = -1 ; x <= 1 ; x++) {
		for (int y = -1 ; y <= 1 ; y++) {
			
			// Current gbuffer depths 
			float SampleDepth = texture(u_PositionTexture, v_TexCoords + vec2(x,y) * TexelSize).x;
			float EffectiveDepth = SampleDepth;
			EffectiveDepth = EffectiveDepth < 0.0f ? 995.0f : EffectiveDepth;
			
			if (EffectiveDepth < MinEffectiveCurrentDepth) {
				MinCurrentDepth = SampleDepth;
				MinEffectiveCurrentDepth = EffectiveDepth;
			}

			if (SampleDepth < 0.0f) { SkyEdge = true; }

			// Previous gbuffer samples
			float SamplePrevDepth = texture(u_PreviousPositionTexture, Reprojected + vec2(x,y) * TexelSize).x;
			float EffectivePrevDepth = SamplePrevDepth;
			EffectivePrevDepth = EffectivePrevDepth < 0.0f ? 995.0f : EffectivePrevDepth;
			
			if (EffectivePrevDepth < MinEffectivePreviousDepth) {
				MinPreviousDepth = SamplePrevDepth;
				MinEffectivePreviousDepth = EffectivePrevDepth;
			}

		}
	}

	return vec2(MinCurrentDepth, MinPreviousDepth);
}

vec4 CatmullRom(sampler2D tex, in vec2 uv)
{
    vec2 texSize = textureSize(tex, 0).xy;
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5f) + 0.5f;
    vec2 f = samplePos - texPos1;
    vec2 w0 = f * (-0.5f + f * (1.0f - 0.5f * f));
    vec2 w1 = 1.0f + f * f * (-2.5f + 1.5f * f);
    vec2 w2 = f * (0.5f + f * (2.0f - 1.5f * f));
    vec2 w3 = f * f * (-0.5f + 0.5f * f);
    
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / (w1 + w2);

    vec2 texPos0 = texPos1 - 1;
    vec2 texPos3 = texPos1 + 2;
    vec2 texPos12 = texPos1 + offset12;

    texPos0 /= texSize;
    texPos3 /= texSize;
    texPos12 /= texSize;

    vec4 result = vec4(0.0f);

    result += texture(tex, vec2(texPos0.x, texPos0.y), 0.0f) * w0.x * w0.y;
    result += texture(tex, vec2(texPos12.x, texPos0.y), 0.0f) * w12.x * w0.y;
    result += texture(tex, vec2(texPos3.x, texPos0.y), 0.0f) * w3.x * w0.y;
    result += texture(tex, vec2(texPos0.x, texPos12.y), 0.0f) * w0.x * w12.y;
    result += texture(tex, vec2(texPos12.x, texPos12.y), 0.0f) * w12.x * w12.y;
    result += texture(tex, vec2(texPos3.x, texPos12.y), 0.0f) * w3.x * w12.y;
    result += texture(tex, vec2(texPos0.x, texPos3.y), 0.0f) * w0.x * w3.y;
    result += texture(tex, vec2(texPos12.x, texPos3.y), 0.0f) * w12.x * w3.y;
    result += texture(tex, vec2(texPos3.x, texPos3.y), 0.0f) * w3.x * w3.y;

    return result;
}

// Samples history texture and applies neighbourhood clipping
vec3 SampleHistory(vec2 Reprojected, vec4 WorldPosition, float Bias) 
{
	const int KernelX = 1;
	const int KernelY = 2;
	
    vec3 MinColor = vec3(100.0);
	vec3 MaxColor = vec3(-100.0); 
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture,0);

	// Find min, max neighbours 
    for(int x = -KernelX; x <= KernelX; x++) {
        for(int y = -KernelY; y <= KernelY; y++) {
            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize).rgb; 
			Sample = Reinhard(Sample);
            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

	MinColor -= Bias;
	MaxColor += Bias;

	vec3 ReprojectedColor = Reinhard(CatmullRom(u_PreviousColorTexture, Reprojected).xyz);
    return (AABBClip((ReprojectedColor), (MinColor), (MaxColor))).xyz;
}

void main()
{
	vec2 Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	vec2 TexCoord = v_TexCoords;
	vec3 CurrentColor = CatmullRom(u_CurrentColorTexture, TexCoord).rgb;

	if (!u_Enabled)
	{
		o_Color = CurrentColor;
		return;
	}

	vec4 WorldPosition = SampleBasePosition(u_PositionTexture).rgba;
	vec3 rD = GetRayDirection(v_TexCoords).xyz;
	rD = normalize(rD);

	float ReduceWeight = 1.0f;
	bool Sky = false;

	// Accumulate fewer frames for the sky because i don't handle cloud motion vectors right now.
	// Temporary!
	if (WorldPosition.w < 0.001f)
	{
		const float rpd = 32.0f;
		WorldPosition.xyz = u_InverseView[3].xyz + rD * rpd;
		WorldPosition.w = rpd;
		ReduceWeight = 0.75f;
		Sky = true;
	}

	vec2 CurrentCoord = v_TexCoords;
	vec2 PreviousCoord;

	// Todo -> Make this use motion vectors 
	vec3 PreviousCameraPos = u_CameraHistory[1];
	vec3 CurrentCameraPos = u_CameraHistory[0];
	vec3 Offset = (PreviousCameraPos - CurrentCameraPos);
	
	bool Moved = length(Offset) > 0.01f;

	bool PannedCamera = u_InversePrevView[0] != u_InverseView[0];
	float PanAmount = distance(u_InversePrevView[0], u_InverseView[0]);
	bool AggressivePan = PanAmount > 0.25f;

	// Reproject 
	PreviousCoord = Reprojection(WorldPosition.xyz);

	float bias = 0.01f;
	
	// Screenspace check
	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		CurrentCoord.x > bias && CurrentCoord.x < 1.0f-bias &&
		CurrentCoord.y > bias && CurrentCoord.y < 1.0f-bias)
	{
		// Tonemap 
		CurrentColor = Reinhard(CurrentColor);

		// History color 

		float ClipBias = 0.06f;

		if (PannedCamera) {
			ClipBias = 0.03f;

			if (AggressivePan) {
				ClipBias = 0.01f;
			}
		}

		if (Moved) {
			ClipBias = 0.0000001f;
		}

		// More aggressive clipping for sky 
		ClipBias *= Sky ? 0.2f : 1.0f;

		vec3 PrevColor = SampleHistory(PreviousCoord, WorldPosition.xyzw, ClipBias);

		// Velocity vector 
		vec2 Velocity = (TexCoord - PreviousCoord.xy) * Dimensions;
		float BlendFactor = exp(-length(Velocity)) * 0.75f + 0.7f;
		BlendFactor = clamp(BlendFactor, 0.0f, 0.975f);

		bool SkySample = false;

		// Depth rejection ->
		if (u_DepthWeight) {
			float BaseExp = 3.0f;

			if (PannedCamera) {
				BaseExp = 6.0f;

				if (AggressivePan) {
					BaseExp = 20.0f;
				}
			}

			if (Moved) {
				BaseExp = 32.0f;
			}

			vec2 MinDepths = GetNearestDepths3x4(PreviousCoord, SkySample);
			bool MovedThresh = distance(u_View[3].xyz, u_PrevView[3].xyz) > 0.001f;
			float DepthRejection = mix(mix(pow(exp(-abs(MinDepths.x - MinDepths.y)), BaseExp * u_DepthExponentMultiplier), 1.0f, (MovedThresh ? 0.0f : 0.5f)), 1.0f, SkySample ? 0.9f : 0.0f);
			DepthRejection = clamp(DepthRejection, 0.0f, 1.0f);
			BlendFactor *= DepthRejection;
		}

		// Avoid clipping artifacts if a block modification is made 
		BlendFactor = u_BlockModified ? max(BlendFactor * 0.5f, 0.225f) : BlendFactor;

		// Blend with history (with reinhard bias)
		o_Color = InverseReinhard(mix((CurrentColor.xyz), (PrevColor.xyz), clamp(BlendFactor * ReduceWeight, 0.001f, 0.95f)));
		
		// Blur factor debug 
		const bool DebugTAA = false;
		if (DebugTAA)
			o_Color = vec3(BlendFactor);


	}

	else 
	{
		o_Color = CurrentColor;
	}
}

