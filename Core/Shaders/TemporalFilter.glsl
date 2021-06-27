#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

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

vec2 Reprojection(vec3 pos) 
{
	return ProjectPositionPrevious(pos).xy * 0.5f + 0.5f;
}

vec4 GetClampedColor(vec2 reprojected, in vec3 worldpos)
{
	int quality = clamp(u_TemporalQuality, 0, 2);

	if (quality == 1)
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

	else if (quality == 2)
	{
		vec2 TexelSize = 1.0f / textureSize(u_PreviousColorTexture, 0).xy;
		vec2 TexelSize2 = 1.0f / textureSize(u_CurrentColorTexture, 0).xy;
		
		vec2 BestOffset = vec2(0.0f, 0.0f);
		float BestDiff = 10000.0f;

		vec4 minclr = vec4(10000.0f);
		vec4 maxclr = vec4(-10000.0f);

		const int BoxSampleSize = 1;

		for(int x = -BoxSampleSize; x <= BoxSampleSize; x++) 
		{
			for(int y = -BoxSampleSize; y <= BoxSampleSize; y++) 
			{
				vec4 SampledPosition = texture(u_PreviousFramePositionTexture, reprojected + (vec2(x, y) * TexelSize)).rgba;
				vec4 Fetch = texture(u_CurrentColorTexture, v_TexCoords + (vec2(x,y) * TexelSize2)).rgba; 

				minclr = min(minclr, Fetch.xyzw); 
				maxclr = max(maxclr, Fetch.xyzw); 

				if (SampledPosition.w > 0.0f)
				{
					float Diff = abs(distance(worldpos, SampledPosition.xyz));

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

		vec4 FinalColor = texture(u_PreviousColorTexture, reprojected + (BestOffset * TexelSize)).xyzw;
		return clamp(FinalColor, minclr, maxclr);
	}

	else 
	{
		return texture(u_PreviousColorTexture, reprojected).rgba;
	}
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;

		if (u_ReflectionTemporal)
		{
			vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
			CameraOffset *= 0.5f;
			Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);
		}

		else 
		{
			Reprojected = Reprojection(CurrentPosition.xyz);
		}

		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec4 PrevColor = GetClampedColor(Reprojected, CurrentPosition.xyz).rgba;
		vec3 PrevPosition = texture(u_PreviousFramePositionTexture, Reprojected).xyz;

		if (Reprojected.x > 0.0 && Reprojected.x < 1.0 && Reprojected.y > 0.0 && Reprojected.y < 1.0)
		{
			float d = abs(distance(PrevPosition, CurrentPosition.xyz));
			float BlendFactor = d;
			BlendFactor = exp(-BlendFactor);
			BlendFactor = clamp(BlendFactor, clamp(u_MinimumMix, 0.05f, 0.9f), clamp(u_MaximumMix, 0.1f, 0.98f));
			
			if (u_ReflectionTemporal)
			{
				//float WeightedDistance = abs(u_CurrentCameraPos.x - u_PrevCameraPos.x) + 
				//					     (abs(u_CurrentCameraPos.y - u_PrevCameraPos.y) * 6.0f) + 
				//						 abs(u_CurrentCameraPos.z - u_PrevCameraPos.z);
				//float BlendMultiplier = 1.0f - (WeightedDistance * 10.0f);
				//BlendMultiplier = clamp(BlendMultiplier, 0.1f, 1.0f);

				const float BlendMultiplier = 0.5f;
				BlendFactor = (u_PrevCameraPos != u_CurrentCameraPos) ? BlendFactor * BlendMultiplier : BlendFactor;
				BlendFactor = clamp(BlendFactor, 0.01f, 0.925f);
				o_Color = mix(CurrentColor, PrevColor, clamp(BlendFactor, 0.0f, u_MaximumMix));
			}

			else 
			{
				o_Color = mix(CurrentColor, PrevColor, BlendFactor);
			}
		}

		else 
		{
			o_Color = mix(CurrentColor, PrevColor, clamp(u_MinimumMix * 0.75f, 0.05f, 0.9f));
		}
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords);
	}
}

