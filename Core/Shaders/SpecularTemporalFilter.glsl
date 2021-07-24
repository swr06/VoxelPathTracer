#version 330 core

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


uniform sampler2D u_PBRTex;

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

vec2 Dimensions;

vec3 ProjectPositionPrevious(vec3 pos)
{
	vec3 WorldPos = pos;
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;

	return ProjectedPosition.xyz;
}

bool Valid(vec2 p)
{
	return p.x > 0.0f && p.x < 1.0f && p.y > 0.0f && p.y < 1.0f;
}

vec2 Reprojection(vec3 pos) 
{
	return ProjectPositionPrevious(pos).xy * 0.5f + 0.5f;
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

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;

		vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
		CameraOffset *= 0.6f;
		Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);

		float RoughnessAt = texture(u_PBRTex, v_TexCoords).r;
		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec3 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyz;
		float d = abs(distance(PrevPosition, CurrentPosition.xyz)); // Disocclusion check
		float ReprojectBias = 0.0125f;

		if (Reprojected.x > 0.0 + ReprojectBias && Reprojected.x < 1.0 - ReprojectBias 
		 && Reprojected.y > 0.0 + ReprojectBias && Reprojected.y < 1.0f - ReprojectBias && d <= 0.5f)
		{
			vec4 PrevColor = texture(u_PreviousColorTexture, Reprojected);
			bool Moved = u_CurrentCameraPos != u_PrevCameraPos;
			float BlendFactor = Moved ? 0.8f : 0.8750f; 
			o_SH = mix(CurrentColor, PrevColor, BlendFactor);

			// store cocg
			vec2 CurrentCoCg = texture(u_CurrentCoCg, v_TexCoords).rg;
			vec2 PrevCocg = texture(u_PrevCoCg, Reprojected).rg;
			o_CoCg = mix(CurrentCoCg, PrevCocg, BlendFactor);
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

