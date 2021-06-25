#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_PositionTexture;

uniform sampler2D u_PreviousColorTexture;

uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

uniform bool u_Enabled;

vec2 View;
vec2 Dimensions;
vec2 TexCoord;

vec2 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

vec3 NeighbourhoodClamping(vec3 color, vec3 tempColor) 
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
		minclr = min(minclr, clr);
		maxclr = max(maxclr, clr);

	}

	return clamp(tempColor, minclr - 0.025f, maxclr + 0.025f);
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;
	TexCoord = v_TexCoords;

	vec3 CurrentColor = texture(u_CurrentColorTexture, TexCoord).rgb;

	if (!u_Enabled)
	{
		o_Color = CurrentColor;
		return;
	}

	vec4 WorldPosition = texture(u_PositionTexture, v_TexCoords).rgba;

	if (WorldPosition.w <= 0.0f)
	{
		o_Color = CurrentColor;
		return;
	}

	vec2 CurrentCoord = v_TexCoords;
	vec2 PreviousCoord = Reprojection(WorldPosition.xyz); 
	float bias = 0.01f;

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		CurrentCoord.x > bias && CurrentCoord.x < 1.0f-bias &&
		CurrentCoord.y > bias && CurrentCoord.y < 1.0f-bias)
	{

		vec3 PrevColor = texture(u_PreviousColorTexture, PreviousCoord).rgb;
		PrevColor.rgb = NeighbourhoodClamping(CurrentColor.rgb, PrevColor.rgb);

		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;

		float BlendFactor = 0.0f;
	
		BlendFactor = exp(-length(velocity)) * 0.2f + 0.6f;
		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, clamp(BlendFactor, 0.025f, 0.9f));
	}

	else 
	{
		o_Color = CurrentColor;
	}
}

