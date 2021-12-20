// Temporal filter/supersampler for specular reflections 

#version 430 core

#define USE_NEW_REPROJECTION 1
#define EPS 0.01f

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture; // -> Current SHy 
uniform sampler2D u_PreviousColorTexture; // -> Previous SHy

uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform sampler2D u_CurrentCoCg;
uniform sampler2D u_PrevCoCg;

uniform sampler2D u_PreviousNormalTexture;


uniform sampler2D u_PBRTex;
uniform sampler2D u_SpecularHitDist;
uniform sampler2D u_PrevSpecularHitDist;

uniform sampler2D u_EmissivityIntersectionMask;

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
uniform bool TEMPORAL_SPEC = false;
uniform bool u_SmartClip = true;

// Firefly rejection ->
uniform bool u_FireflyRejection;
uniform bool u_AggressiveFireflyRejection;

vec2 Dimensions;

bool Valid(vec2 p)
{
	return p.x > 0.0f && p.x < 1.0f && p.y > 0.0f && p.y < 1.0f;
}

vec2 Reprojection(vec3 pos) 
{
	vec3 WorldPos = pos;
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xy * 0.5f + 0.5f;
}

vec2 Project(vec3 pos) 
{
	vec3 WorldPos = pos;
	vec4 ProjectedPosition = u_Projection * u_View * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xy * 0.5f + 0.5f;
}

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
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

vec3 SampleNormal(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}

vec3 SHToIrridiance(vec4 shY, vec2 CoCg)
{
    float Y = max(0, 3.544905f * shY.w);
    Y = max(Y, 0.0);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

float remap(float original_value, float original_min, float original_max, float new_min, float new_max)
{
    return new_min + (clamp((original_value - original_min) / (original_max - original_min), 0.0f, 1.0f) * (new_max - new_min));
}

vec3 ClipToAABB(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

void ReflectionClipping(inout vec4 PreviousSH, inout vec2 PreviousCoCg, vec2 Reprojected, float Roughness) {
	
	const float RoughnessThreshold = 0.275f + 0.01f;
	
	if (Roughness > 0.5f + 0.01f) {
		return;
	}
	
	vec4 MinSH = vec4(1000.0f), MaxSH = vec4(-1000.0f);
	vec2 MinCoCg = vec2(1000.0f), MaxCoCg = vec2(-1000.0f);
	vec3 MinRadiance = vec3(1000.0f);
	vec3 MaxRadiance = vec3(-1000.0f);

	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	
	float AdditionalMaxBias = 0.0f;

	for (int x = -1; x <= 1 ; x++) {
		for (int y = -1 ; y <= 1 ; y++) {
			
			vec2 SampleCoord = Reprojected + vec2(x,y) * TexelSize;
			float Mask = texture(u_EmissivityIntersectionMask, SampleCoord).x;

			bool MaskE = Mask < 0.01f;
			
			if (!MaskE) {
				AdditionalMaxBias += 0.1f;
			}

			vec4 SampleSH = texture(u_CurrentColorTexture, SampleCoord);
			vec2 SampleCoCg = texture(u_CurrentCoCg, SampleCoord).xy;
			
			vec3 Irradiance = SHToIrridiance(SampleSH, SampleCoCg);
			
			if (MinRadiance != min(Irradiance, MinRadiance)) {
				MinSH = SampleSH;
				MinCoCg = SampleCoCg;
				MinRadiance = min(Irradiance, MinRadiance);
			}
			
			if (MaxRadiance != max(Irradiance, MaxRadiance)) {
				MaxSH = SampleSH;
				MaxCoCg = SampleCoCg;
				MaxRadiance = max(Irradiance, MaxRadiance);
			}
		}
	}

	bool ADD_BIAS = true;
	
	if (ADD_BIAS) {
	
		
		bool Smoothish = Roughness < RoughnessThreshold;
		bool Roughish = Roughness > RoughnessThreshold && Roughness < 0.50f + 0.01f;
		
		if (Smoothish) {
		
			float PercievedRoughness = Roughness * Roughness;
			const float RoughnessThresholdSquared = RoughnessThreshold * RoughnessThreshold;
			
			float RemappedRoughness = remap(PercievedRoughness, 0.0f, RoughnessThresholdSquared, 0.0f, 1.0f);
			
			float Bias = mix(0.01f, 0.085f, RemappedRoughness);
			
			if (Roughness > 0.1925f) {
				Bias *= 1.75f;
			}
			
			MinRadiance -= (Bias * 0.95f);
			MaxRadiance += (Bias * 0.95f) + AdditionalMaxBias;
		}
		
		else if (Roughish) 
		{
			MinRadiance -= 0.37f;
			MaxRadiance += 0.37f + (AdditionalMaxBias * 1.1f);
		}
	}

	// Clipping ->
	vec3 BaseRadiance = SHToIrridiance(PreviousSH, PreviousCoCg);
	vec3 ClampedRadiance = ClipToAABB(BaseRadiance, MinRadiance, MaxRadiance);
	//vec3 ClampedRadiance = clamp(BaseRadiance, MinRadiance, MaxRadiance);

	if (ClampedRadiance != BaseRadiance) {
		if (distance(ClampedRadiance, MinRadiance) > distance(ClampedRadiance, MaxRadiance)) {
			PreviousSH = MaxSH;
			PreviousCoCg = MaxCoCg;
		}

		else {
			PreviousSH = MinSH;
			PreviousCoCg = MinCoCg; 
		}
	}
	
	
	return;
}

float SHToY(vec4 shY)
{
    return max(0, 3.544905f * shY.w); // get luminance (Y) from the first spherical harmonic band <-
}

float SHToY(float shY)
{
    return max(0, 3.544905f * shY); // get luminance (Y) from the first spherical harmonic band <-
}


vec3 VarianceFireflyRejection(vec3 Radiance, float VarianceEstimate, vec3 Mean)
{
    vec3 StandardDeviation = vec3(sqrt(max(0.00001f, VarianceEstimate))); // Calculate standard deviation 
    vec3 Threshold = 0.1f + Mean + StandardDeviation * 8.0;
    vec3 ErrorEstimate = vec3(max(vec3(0.0f), Radiance - Threshold));
    return clamp(Radiance - ErrorEstimate, 0.0f, 16.0f); 
}

void FireflyReject(inout vec4 InputSH, inout vec2 InputCoCg) 
{
	float Y = SHToY(InputSH);
	int SampleThreshold = u_AggressiveFireflyRejection ? 3 : 4;

	vec2 Offsets[4] = vec2[4](vec2(1.0f, 0.0f), vec2(0.0f, 1.0f), vec2(-1.0f, 0.0f), vec2(0.0f, -1.0f));
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);

	vec4 NonLitSHySum = vec4(0.0f);
	vec2 NonLitCoCgSum = vec2(0.0f);
	float MinY = 1000.0f;
	int UnlitSampleSum = 0;
	
	for (int i = 0 ; i < 4 ; i++) {
		vec2 SampleCoord = v_TexCoords + Offsets[i] * TexelSize;
		float Mask = texture(u_EmissivityIntersectionMask, SampleCoord).x;
		vec4 SHy = texture(u_CurrentColorTexture, SampleCoord);

		if (Mask < 0.01f) {
			NonLitSHySum += SHy;
			NonLitCoCgSum += texture(u_CurrentCoCg, SampleCoord).xy;
			UnlitSampleSum++;
		}

		float yy = SHToY(SHy);
		MinY = min(MinY, yy);
	}

	if (UnlitSampleSum >= SampleThreshold) {

		// Average ->
		NonLitSHySum /= float(UnlitSampleSum);
		NonLitCoCgSum /= float(UnlitSampleSum);

		InputSH = NonLitSHySum;
		InputCoCg = NonLitCoCgSum;
	}

	return;
}

void FireflyRejectLuminosity(inout vec4 InputSH, inout vec2 InputCoCg) 
{
	float Y = SHToY(InputSH);

	vec2 Offsets[4] = vec2[4](vec2(1.0f, 0.0f), vec2(0.0f, 1.0f), vec2(-1.0f, 0.0f), vec2(0.0f, -1.0f));
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	int SampleThreshold = u_AggressiveFireflyRejection ? 3 : 4;

	vec4 NonLitSHy = vec4(0.0f);
	vec2 NonLitCoCg = vec2(0.0f);
	float MinY = 1000.0f;
	int Count = 0;
	
	for (int i = 0 ; i < 4 ; i++) {
		vec2 SampleCoord = v_TexCoords + Offsets[i] * TexelSize;
		vec4 SHy = texture(u_CurrentColorTexture, SampleCoord);

		float yy = SHToY(SHy);
		float yerror = abs(yy - Y);

		if (yerror > 0.25f) {
			NonLitSHy = SHy;
			NonLitCoCg = texture(u_CurrentCoCg, SampleCoord).xy;
			Count++;
		}

		MinY = min(MinY, yy);
	}

	if (Count >= SampleThreshold) {
		InputSH = NonLitSHy;
		InputCoCg = NonLitCoCg;
	}

	return;
}



void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;
	vec3 InitialNormal = SampleNormal(u_NormalTexture, v_TexCoords);
	vec4 CurrentSHy = texture(u_CurrentColorTexture, CurrentCoord).rgba;
	vec2 CurrentCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;
	
	// Firefly rejection ->

	if (u_FireflyRejection) {
		float Mask = texture(u_EmissivityIntersectionMask, v_TexCoords).x;

		if (Mask > 0.05f) 
		{
			FireflyReject(CurrentSHy, CurrentCoCg);
		}
	}

	if (CurrentPosition.a > 0.0f && TEMPORAL_SPEC)
	{
		float HitDistanceCurrent = texture(u_SpecularHitDist, v_TexCoords).r;
		bool SkySample = HitDistanceCurrent < 0.0f;

		float RoughnessAt = texture(u_PBRTex, v_TexCoords).r;
		
		bool LessValid = false;
		vec2 Reprojected; 

		const bool UseNewReprojection = bool(USE_NEW_REPROJECTION);
		
		// Usually we wouldn't be able to use such a high roughness threshold 
		// But this works since we apply a bilateral filter to the distance reproj data
		
		
		
		if (RoughnessAt <= 0.70f+0.01f && HitDistanceCurrent > 0.0f && UseNewReprojection && !SkySample)
		{
			// Reconstruct the reflected position to properly reproject
			vec3 I = normalize(v_RayOrigin - CurrentPosition.xyz);
			vec3 ReflectedPosition = CurrentPosition.xyz - I * HitDistanceCurrent;

			// Project and use the approximated motion vector :
			vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(ReflectedPosition, 1.0f);
			ProjectedPosition.xyz /= ProjectedPosition.w;
			Reprojected = ProjectedPosition.xy * 0.5f + 0.5f;

			// Validate hit distance : 

			float PreviousT = texture(u_PrevSpecularHitDist, Reprojected).x;
			
			if (abs(PreviousT - HitDistanceCurrent) >= 2.2f) {
			
				// If hit distance is invalid, fallback on very-approximate reprojection
				
				vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
				CameraOffset *= 0.6f;
				Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);
				LessValid = true;
			}
		}

		else {
			// very approximate reprojection based on camera position delta
			
			vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
			CameraOffset *= 0.6f;
			Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);
			LessValid = true;
		}

		if (SkySample) {
			
			// We can make some reprojection assumptions when the sky is reflected (such as a constant transversal) ->
			float SkyTransversalApproximate = 64.0f;
			vec3 I = normalize(v_RayOrigin - CurrentPosition.xyz);
			vec3 ReflectedPosition = CurrentPosition.xyz - I * SkyTransversalApproximate;
			vec4 ProjectedPosition = u_PrevProjection * u_PrevView  * vec4(ReflectedPosition, 1.0f);
			ProjectedPosition.xyz /= ProjectedPosition.w;
			Reprojected = ProjectedPosition.xy * 0.5f + 0.5f;
		}
		
		// Surface disocclusion check : 
		vec4 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyzw;
		float d = abs(distance(PrevPosition.xyz, CurrentPosition.xyz)); // Disocclusion check

		float ReprojectBias = 0.01f;
		
		vec3 PrevNormal = SampleNormal(u_PreviousNormalTexture, Reprojected.xy);

		if (Reprojected.x > 0.0 + ReprojectBias && Reprojected.x < 1.0 - ReprojectBias 
		 && Reprojected.y > 0.0 + ReprojectBias && Reprojected.y < 1.0f - ReprojectBias && 
		 d < 0.85f && PrevNormal==InitialNormal)
		{
			vec4 PrevSH;
			vec2 PrevCoCg;
			//ComputeClamped(Reprojected, PrevSH, PrevCoCg);

			PrevSH = texture(u_PreviousColorTexture, Reprojected);
			PrevCoCg = texture(u_PrevCoCg, Reprojected).xy;
			
			vec4 PrevSHBck = PrevSH;
			vec2 PrevCoCgBck = PrevCoCg;

			const float DistanceEps = 0.0025f;
			bool Moved = distance(u_CurrentCameraPos, u_PrevCameraPos) > DistanceEps;
			//bool Moved = u_CurrentCameraPos != u_PrevCameraPos;
			
			// Compute accumulation factor ->
			float BlendFactor = LessValid ? (Moved ? 0.725f : 0.875f) : (Moved ? 0.875f : 0.95f); 

			//bool FuckingSmooth = RoughnessAt <= 0.25f + 0.01f;
			bool TryClipping = RoughnessAt < 0.5f + 0.01f;

			// Radiance clipping -> 
			// Used to reduce ghosting 
			if (TryClipping && Moved && u_SmartClip) {

				// Clip sample ->
				ReflectionClipping(PrevSH, PrevCoCg, v_TexCoords, RoughnessAt);
				
				
				if (RoughnessAt > 0.26f) 
				{
					PrevSH = mix(PrevSHBck, PrevSH, 0.9f);
					PrevCoCg = mix(PrevCoCgBck, PrevCoCg, 0.9f);
					BlendFactor *= 1.072525f;
				}
				
			
			}

			bool PreviousSkySample = texture(u_PrevSpecularHitDist, Reprojected.xy).x < 0.0f;

			//if (SkySample == PreviousSkySample) {
			//	BlendFactor *= 1.15f;
			//}

			BlendFactor = clamp(BlendFactor, 0.001f, 0.95f);

			// mix sh
			o_SH = mix(CurrentSHy, PrevSH, BlendFactor);
			o_CoCg = mix(CurrentCoCg, PrevCoCg, BlendFactor);
		}

		else 
		{
			o_SH = CurrentSHy.xyzw;
			o_CoCg = CurrentCoCg;
		}
	}

	else 
	{
		o_SH = CurrentSHy;
		o_CoCg = CurrentCoCg;
	}
}

