#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

uniform float u_MixModifier = 0.8;

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

// Reduce smearing 
void GetNearestDepth(in float input_depth, out float odepth, out vec2 best_offset, vec2 reprojected) 
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

	odepth = 100000.0f;

	for(int i = 0; i < 8; i++) 
	{
		vec2 offset = neighbourhoodOffsets[i] * View;

		float depth_at = texture(u_PreviousFramePositionTexture, reprojected + offset).r;
		float d = abs(input_depth - depth_at);

		if (d < odepth)
		{
			odepth = d;
			best_offset = neighbourhoodOffsets[i]; 
		}
	}
}


vec4 GetClampedColor(vec2 reprojected)
{
	ivec2 Coord = ivec2(v_TexCoords * Dimensions); 

	vec4 minclr = vec4(10000.0f); 
	vec4 maxclr = vec4(-10000.0f); 

	for(int x = -2; x <= 2; x++) 
	{
		for(int y = -2; y <= 2; y++) 
		{
			vec4 Fetch = texelFetch(u_CurrentColorTexture, Coord + ivec2(x,y), 0); 

			minclr = min(minclr, Fetch); 
			maxclr = max(maxclr, Fetch); 
		}
	}

	minclr -= 0.075f; 
	maxclr += 0.075f; 
	
	return clamp(texture(u_PreviousColorTexture, reprojected), minclr, maxclr); 

}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;

	TexCoord = v_TexCoords;

	vec2 CurrentCoord = TexCoord;
	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;
	vec3 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgb;

	if (CurrentPosition.a > 0.0f && CurrentColor.x > -0.9f && CurrentColor.y > -0.9 && CurrentColor.z > -0.9)
	{
		vec2 PreviousCoord = Reprojection(CurrentPosition.xyz); 

		vec3 PrevColor = GetClampedColor(PreviousCoord).rgb;

		vec3 AverageColor;
		float ClosestDepth;

		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;

		float BlendFactor = float(
			PreviousCoord.x > 0.0 && PreviousCoord.x < 1.0 &&
			PreviousCoord.y > 0.0 && PreviousCoord.y < 1.0
		);

		//float CurrentDepth = CurrentPosition.z;
		//float PreviousDepth = texture(u_PreviousFramePositionTexture, PreviousCoord).z;
		//float DepthDifference = abs(CurrentDepth - PreviousDepth);

		BlendFactor *= exp(-length(velocity)) * 0.35f;
		BlendFactor += u_MixModifier;
		BlendFactor = clamp(BlendFactor, 0.01f, 0.9790f);
		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, BlendFactor);
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgb;
	}
}

