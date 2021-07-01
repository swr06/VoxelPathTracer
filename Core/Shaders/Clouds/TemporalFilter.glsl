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

uniform vec3 u_CurrentPosition;
uniform vec3 u_PreviousPosition;

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

float FastDist(vec3 p1, vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

float GetLuminance(vec3 color)
{
	return dot(color, vec3(0.299f, 0.587f, 0.114f));
}

vec4 GetBestSample(vec2 reprojected, vec3 CurrPos)
{
	vec2 TexelSize = 1.0f / textureSize(u_PreviousColorTexture, 0).xy;
	vec2 BestOffset = vec2(0.0f, 0.0f);
	float BestDiff = 10000.0f;
	vec4 minclr = vec4(10000.0f);
	vec4 maxclr = vec4(-10000.0f);
	const int BoxSampleSize = 2; 

	// 4*4 = 16 sample box
	for(int x = -BoxSampleSize; x <= BoxSampleSize; x++) 
	{
		for(int y = -BoxSampleSize; y <= BoxSampleSize; y++) 
		{
			vec4 SampledPosition = texture(u_PreviousFramePositionTexture, reprojected + (vec2(x, y) * TexelSize)).rgba;
			vec4 Fetch = texture(u_CurrentColorTexture, v_TexCoords + (vec2(x,y) * TexelSize)).rgba; 

			minclr = min(minclr, Fetch.xyzw); 
			maxclr = max(maxclr, Fetch.xyzw); 

			if (SampledPosition.w > 0.0f)
			{
				float Diff = abs(FastDist(CurrPos, SampledPosition.xyz));

				if (Diff < BestDiff)
				{
					BestDiff = Diff;
					BestOffset = vec2(x, y);
				}
			}	
		}
	}

	minclr -= 0.010f; 
	maxclr += 0.010f; 

	vec4 FinalColor = texture(u_PreviousColorTexture, reprojected + (BestOffset * TexelSize)).xyzw;
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
		vec3 CameraOffset = u_CurrentPosition - u_PreviousPosition;
		CameraOffset *= 0.0;
		vec2 PreviousCoord = Reprojection(CurrentPosition.xyz + CameraOffset); 
		vec4 PrevColor = GetBestSample(PreviousCoord, CurrentPosition.xyz);
		vec4 CurrentColor = texture(u_CurrentColorTexture, v_TexCoords).rgba;
		vec4 PrevPosition = texture(u_PreviousFramePositionTexture, PreviousCoord).xyzw;

		vec3 AverageColor;
		float ClosestDepth;
		float error = distance(PrevPosition.xyz, CurrentPosition.xyz);
		const float error_thresh = 0.41414f;

		if(PreviousCoord.x > 0.0 && PreviousCoord.x < 1.0 && PreviousCoord.y > 0.0 && PreviousCoord.y < 1.0 && 
		   error < error_thresh && PrevPosition.a != 0.0f)
		{
			o_Color = mix(CurrentColor.xyzw, PrevColor.xyzw, 0.95f);
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
