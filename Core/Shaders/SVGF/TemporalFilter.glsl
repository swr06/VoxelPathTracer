#version 330 core

layout (location = 0) out vec4 o_Color;
layout (location = 1) out vec2 o_SH;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform sampler2D u_PreviousSH;
uniform sampler2D u_CurrentSH;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_MinimumMix = 0.25f;
uniform float u_MaximumMix = 0.975f;
uniform vec3 u_PrevCameraPos;
uniform vec3 u_CurrentCameraPos;

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

void main()
{
	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;
		Reprojected = Reprojection(CurrentPosition.xyz);
		
		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec4 PrevColor = texture(u_PreviousColorTexture, Reprojected);
		vec3 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyz;

		float Bias = 0.006524f;

		if (Reprojected.x > 0.0 + Bias && Reprojected.x < 1.0 - Bias && Reprojected.y > 0.0 + Bias && Reprojected.y < 1.0 - Bias)
		{
			float d = abs(distance(PrevPosition, CurrentPosition.xyz));
			float BlendFactor = d;
			BlendFactor = exp(-BlendFactor);
			BlendFactor = clamp(BlendFactor, clamp(u_MinimumMix, 0.01f, 0.9f), clamp(u_MaximumMix, 0.1f, 0.98f));
			o_Color = mix(CurrentColor, PrevColor, BlendFactor);
			vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
			vec2 PrevSH = texture(u_PreviousSH, Reprojected.xy).xy;
			o_SH = mix(CurrentSH, PrevSH, BlendFactor);
		}

		else 
		{
			o_Color = CurrentColor;
			vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
			o_SH = CurrentSH;
		}
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords);
		vec2 CurrentSH = texture(u_CurrentSH, v_TexCoords).xy;
		o_SH = CurrentSH;
	}
}

