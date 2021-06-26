#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_PositionTexture;

uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousPositionTexture;

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

float FastDistance(in vec3 p1, in vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

vec2 GetNearFragment(vec2 coord, vec3 posat)
{
	vec2 NearestOffset = vec2(0.0f, 0.0f);
	float Difference = 10000;

	vec2 Offsets[9] = vec2[9]
	(
		vec2(-1.0, -1.0),
		vec2( 0.0, -1.0),
		vec2( 1.0, -1.0),
		vec2(-1.0,  0.0),
		vec2(0.0f, 0.0f),
		vec2( 1.0,  0.0),
		vec2(-1.0,  1.0),
		vec2( 0.0,  1.0),
		vec2( 1.0,  1.0)
	);

	vec2 TexelSize = 1.0f / textureSize(u_PositionTexture, 0);
	vec2 TexelSize2 = 1.0f / textureSize(u_PreviousColorTexture, 0);
	
	for (int i = 0 ; i < 9 ; i++)
	{
		vec4 PositionAt = texture(u_PreviousPositionTexture, coord + (Offsets[i] * TexelSize));

		if (PositionAt.w > 0.0f) 
		{
			float diff = abs(FastDistance(posat, PositionAt.xyz));

			if (diff < Difference)
			{
				NearestOffset = Offsets[i];
				Difference = diff;
			}
		}
	}

	return coord + NearestOffset * TexelSize2;
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
	PreviousCoord = GetNearFragment(PreviousCoord, WorldPosition.xyz);

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		CurrentCoord.x > bias && CurrentCoord.x < 1.0f-bias &&
		CurrentCoord.y > bias && CurrentCoord.y < 1.0f-bias)
	{
		vec4 WorldPositionPrev = texture(u_PositionTexture, PreviousCoord).rgba;

		if (WorldPositionPrev.w <= 0.0f)
		{
			o_Color = CurrentColor;
			return;
		}

		vec3 PrevColor = texture(u_PreviousColorTexture, PreviousCoord).rgb;
		PrevColor = NeighbourhoodClamping(CurrentColor, PrevColor);

		// Construct our motion vector
		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;
		float BlendFactor = exp(-length(velocity)) * 0.2f + 0.6f;
		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, clamp(BlendFactor, 0.025f, 0.925f));
	}

	else 
	{
		o_Color = CurrentColor;
	}
}

