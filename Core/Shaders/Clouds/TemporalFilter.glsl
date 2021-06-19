#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

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

float GetLuminance(vec3 color)
{
	return dot(color, vec3(0.299f, 0.587f, 0.114f));
}

vec3 GetBestSample(vec2 reprojected, vec3 world_pos)
{
	vec2 TexelSize = 1.0f / textureSize(u_PreviousColorTexture, 0).xy;
	vec2 TexelSize2 = 1.0f / textureSize(u_CurrentColorTexture, 0).xy;

	vec2 BestOffset = vec2(0.0f, 0.0f);
	float BestDiff = 10000.0f;

	vec3 minclr = vec3(10000.0f);
	vec3 maxclr = vec3(-10000.0f);

	const int BoxSampleSize = 1;

	for(int x = -BoxSampleSize; x <= BoxSampleSize; x++) 
	{
		for(int y = -BoxSampleSize; y <= BoxSampleSize; y++) 
		{
			vec4 SampledPosition = texture(u_PreviousFramePositionTexture, reprojected + (vec2(x, y) * TexelSize)).rgba;
			vec3 Fetch = texture(u_CurrentColorTexture, v_TexCoords + (vec2(x,y) * TexelSize2)).rgb; 

			minclr = min(minclr, Fetch.xyz); 
			maxclr = max(maxclr, Fetch.xyz); 

			if (SampledPosition.w > 0.0f)
			{
				float Diff = abs(distance(world_pos, SampledPosition.xyz));

				if (Diff < BestDiff)
				{
					BestDiff = Diff;
					BestOffset = vec2(x, y);
				}
			}	
		}
	}

	minclr -= 0.065f; 
	maxclr += 0.065f; 

	vec3 FinalColor = texture(u_PreviousColorTexture, reprojected + (BestOffset * TexelSize)).xyz;
	return clamp(FinalColor, minclr, maxclr);
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;

	TexCoord = v_TexCoords;

	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a != 0.0f)
	{
		vec2 PreviousCoord = Reprojection(CurrentPosition.xyz); 
		vec3 PrevColor = GetBestSample(PreviousCoord, CurrentPosition.xyz).rgb;
		vec3 CurrentColor = texture(u_CurrentColorTexture, v_TexCoords).rgb;

		vec3 AverageColor;
		float ClosestDepth;

		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;

		if(PreviousCoord.x > 0.0 && PreviousCoord.x < 1.0 && PreviousCoord.y > 0.0 && PreviousCoord.y < 1.0)
		{
			float BlendFactor = 1.0f;
			BlendFactor = exp(-length(velocity)) * 0.2f;
			BlendFactor += 0.8f;
			BlendFactor = clamp(BlendFactor, 0.925f, 0.98f); 
			o_Color = mix(CurrentColor.xyz, PrevColor.xyz, BlendFactor);
		}

		else 
		{
			o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgb;
		}
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords).rgb;
	}
}
