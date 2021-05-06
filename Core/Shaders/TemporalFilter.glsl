#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

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

vec3 NeighbourhoodClamping(vec3 color, vec3 tempColor, out vec3 avg_col) 
{
	vec2 neighbourhoodOffsets[8] = vec2[8]
	(
		vec2(-1.0, -1.0),
		vec2( 0.0, -1.0),
		vec2( 1.0, -1.0),
		vec2(-1.0,  0.0),
		vec2( 1.0,  0.0),
		vec2(-1.0,  1.0),
		vec2( 0.0,  1.0),
		vec2( 1.0,  1.0)
	);

	vec3 minclr = color, maxclr = color;

	for(int i = 0; i < 8; i++) 
	{
		vec2 offset = neighbourhoodOffsets[i] * View;
		vec3 clr = texture(u_CurrentColorTexture, TexCoord + offset, 0.0).rgb;
		avg_col += clr;
		minclr = min(minclr, clr);
		maxclr = max(maxclr, clr);
	}

	avg_col /= 8.0f;

	return clamp(tempColor, minclr, maxclr);
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;

	TexCoord = v_TexCoords;

	vec2 CurrentCoord = TexCoord;
	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 PreviousCoord = Reprojection(CurrentPosition.xyz); 

		vec3 CurrentColor = texture(u_CurrentColorTexture, TexCoord).rgb;
		vec3 PrevColor = texture(u_PreviousColorTexture, PreviousCoord).rgb;

		vec3 AverageColor;
		float ClosestDepth;
		vec2 BestOffset;

		//PrevColor.rgb = NeighbourhoodClamping(CurrentColor.rgb, PrevColor.rgb, AverageColor);

		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;

		float BlendFactor = float(
			PreviousCoord.x > 0.0 && PreviousCoord.x < 1.0 &&
			PreviousCoord.y > 0.0 && PreviousCoord.y < 1.0
		);

		BlendFactor *= (exp(-length(velocity)) * 0.62f) + 0.55f; // 0.35f
		BlendFactor = clamp(BlendFactor, 0.03, 0.95);

		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, BlendFactor);
		//o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgb;
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgb;
	}
}

