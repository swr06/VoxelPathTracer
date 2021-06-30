#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

uniform float u_MixModifier = 0.8;
uniform float u_Time;

vec2 View;
vec2 Dimensions;
vec2 TexCoord;

vec2 Reprojection(vec3 pos) 
{
	vec3 WorldPos = pos;

	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

float GetLuminance(vec3 color)
{
	return dot(color, vec3(0.299f, 0.587f, 0.114f));
}


void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;

	TexCoord = v_TexCoords;

	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a != 0.0f)
	{
		vec3 time = vec3(u_Time, 0.0f, u_Time * 0.5f);
		time *= 0.00250f; 
		time *= 0.0f;

		vec2 PreviousCoord = Reprojection(CurrentPosition.xyz + time); 
		vec4 PrevColor = texture(u_PreviousColorTexture, PreviousCoord).rgba;
		vec4 CurrentColor = texture(u_CurrentColorTexture, v_TexCoords).rgba;

		vec3 AverageColor;
		float ClosestDepth;

		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;

		if(PreviousCoord.x > 0.0 && PreviousCoord.x < 1.0 && PreviousCoord.y > 0.0 && PreviousCoord.y < 1.0)
		{
			float BlendFactor = 1.0f;
			BlendFactor = exp(-length(velocity)) * 0.2f;
			BlendFactor += 0.8f;
			BlendFactor = 0.9; 
			o_Color = mix(CurrentColor.xyzw, PrevColor.xyzw, 0.9);
		}

		else 
		{
			o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgba;
		}
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgba;
	}
}
