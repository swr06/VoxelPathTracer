#version 330 core

#define USE_NEW_REPROJECTION 1

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform sampler2D u_CurrentCoCg;
uniform sampler2D u_PrevCoCg;

uniform sampler2D u_PreviousNormalTexture;


uniform sampler2D u_PBRTex;
uniform sampler2D u_SpecularHitDist;
uniform sampler2D u_PrevSpecularHitDist;

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

void ComputeClamped(vec2 r, out vec4 sh, out vec2 cocg) {
	
	vec4 MinSH = vec4(1000.0f), MaxSH = vec4(-1000.0f);
	vec2 MinCoCg = vec2(1000.0f), MaxCoCg = vec2(-1000.0f);

	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);

	for (int x = -1; x <= 1 ; x++) {
		for (int y = -1 ; y <= 1 ; y++) {
			
			vec2 SampleCoord = r + vec2(x,y) * TexelSize;
			vec4 SampleSH = texture(u_CurrentColorTexture, SampleCoord);
			vec2 SampleCoCg = texture(u_CurrentCoCg, SampleCoord).xy;
			MinSH = min(SampleSH, MinSH);
			MaxSH = max(SampleSH, MaxSH);
			MinCoCg = min(SampleCoCg, MinCoCg);
			MaxCoCg = max(SampleCoCg, MaxCoCg);
		}
	}

	float bias = 0.0f;
	MinCoCg -= bias; MinSH -= bias;
	MaxCoCg += bias; MaxSH += bias;

	sh = clamp(texture(u_PreviousColorTexture, r), MinSH, MaxSH);
	cocg = clamp(texture(u_PrevCoCg, r).xy, MinCoCg, MaxCoCg);
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;
	vec3 InitialNormal = SampleNormal(u_NormalTexture, v_TexCoords);

	if (CurrentPosition.a > 0.0f&&TEMPORAL_SPEC)
	{
		float HitDistanceCurrent = texture(u_SpecularHitDist, v_TexCoords).r;

		float RoughnessAt = texture(u_PBRTex, v_TexCoords).r;
		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;

		bool LessValid = false;
		vec2 Reprojected; 

		const bool UseNewReprojection = bool(USE_NEW_REPROJECTION);
		
		if (RoughnessAt < 0.35f && HitDistanceCurrent > 0.0f && UseNewReprojection)
		{
			// Reconstruct the reflected position to properly reproject
			vec3 I = normalize(v_RayOrigin - CurrentPosition.xyz);
			vec3 ReflectedPosition = CurrentPosition.xyz - I * HitDistanceCurrent;

			// Project and use the motion vector :
			vec4 ProjectedPosition = u_PrevProjection * u_PrevView  * vec4(ReflectedPosition, 1.0f);
			ProjectedPosition.xyz /= ProjectedPosition.w;
			Reprojected = ProjectedPosition.xy * 0.5f + 0.5f;

			// Validate hit distance : 

			float PreviousT = texture(u_PrevSpecularHitDist, Reprojected).x;
			if (abs(PreviousT - HitDistanceCurrent) > 0.4f) {
				vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
				CameraOffset *= 0.6f;
				Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);
				LessValid = true;
			}
		}

		else {
			vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
			CameraOffset *= 0.6f;
			Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);
			LessValid = true;
		}
		
		// Disocclusion check : 
		vec3 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyz;
		float d = abs(distance(PrevPosition, CurrentPosition.xyz)); // Disocclusion check

		float ReprojectBias = 0.01f;

		vec3 PrevNormal = SampleNormal(u_PreviousNormalTexture, Reprojected.xy);

		if (Reprojected.x > 0.0 + ReprojectBias && Reprojected.x < 1.0 - ReprojectBias 
		 && Reprojected.y > 0.0 + ReprojectBias && Reprojected.y < 1.0f - ReprojectBias && 
		 d < 0.75f && PrevNormal==InitialNormal)
		{
			vec4 PrevSH;
			vec2 PrevCoCg;
			//ComputeClamped(Reprojected, PrevSH, PrevCoCg);

			PrevSH = texture(u_PreviousColorTexture, Reprojected);
			PrevCoCg = texture(u_PrevCoCg, Reprojected).xy;

			bool Moved = u_CurrentCameraPos != u_PrevCameraPos;
			float BlendFactor = LessValid ? (Moved ? 0.725f : 0.85f) : (Moved ? 0.8f : 0.9f); 

			// mix sh
			o_SH = mix(CurrentColor, PrevSH, BlendFactor);

			// store cocg
			vec2 CurrentCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;
			o_CoCg = mix(CurrentCoCg, PrevCoCg, BlendFactor);
		}

		else 
		{
			o_SH = CurrentColor.xyzw;
			vec2 CurrentCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;
			o_CoCg = CurrentCoCg;
		}
	}

	else 
	{
		o_SH = texture(u_CurrentColorTexture, v_TexCoords);
		vec2 CurrentCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;
		o_CoCg = CurrentCoCg;
	}
}

