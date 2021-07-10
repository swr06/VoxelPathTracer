#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

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
uniform float u_ClampBias = 0.025f;

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

vec4 GetClampedColor(vec2 reprojected, in vec3 worldpos)
{
	int quality = clamp(u_TemporalQuality, 0, 2);
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	ivec2 Coord = ivec2(v_TexCoords * Dimensions); 

	if (quality == 1)
	{

		vec4 minclr = vec4(10000.0f); 
		vec4 maxclr = vec4(-10000.0f); 

		const int SampleSize = 1;

		for(int x = -SampleSize; x <= SampleSize; x++) 
		{
			for(int y = -SampleSize; y <= SampleSize; y++) 
			{
				vec4 Fetch = texelFetch(u_CurrentColorTexture, Coord + ivec2(x,y), 0); 
				minclr = min(minclr, Fetch); 
				maxclr = max(maxclr, Fetch); 
			}
		}

		minclr -= u_ClampBias; 
		maxclr += u_ClampBias; 
		
		return clamp(texture(u_PreviousColorTexture, reprojected), minclr, maxclr); 
	}

	else if (quality == 2)
	{
		vec2 TexelSize = 1.0f / textureSize(u_PreviousColorTexture, 0).xy;
		
		vec2 BestOffset = vec2(0.0f, 0.0f);
		float BestDiff = 10000.0f;

		vec4 minclr = vec4(10000.0f);
		vec4 maxclr = vec4(-10000.0f);

		const int BoxSampleSize = 2;

		for(int x = -BoxSampleSize; x <= BoxSampleSize; x++) 
		{
			for(int y = -BoxSampleSize; y <= BoxSampleSize; y++) 
			{
				vec4 SampledPosition = GetPositionAt(u_PreviousFramePositionTexture, reprojected + (vec2(x, y) * TexelSize)).rgba;
				vec4 Fetch = texelFetch(u_CurrentColorTexture, Coord + ivec2(x,y), 0); 

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

		minclr -= u_ClampBias; 
		maxclr += u_ClampBias; 

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
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;
		Reprojected = Reprojection(CurrentPosition.xyz);

		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec4 PrevColor = GetClampedColor(Reprojected, CurrentPosition.xyz).rgba;
		vec3 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyz;

		if (Reprojected.x > 0.0 && Reprojected.x < 1.0 && Reprojected.y > 0.0 && Reprojected.y < 1.0)
		{
			float d = abs(distance(PrevPosition, CurrentPosition.xyz));
			float BlendFactor = d;
			BlendFactor = exp(-BlendFactor);
			BlendFactor = clamp(BlendFactor, clamp(u_MinimumMix, 0.01f, 0.9f), clamp(u_MaximumMix, 0.1f, 0.98f));
			o_Color = mix(CurrentColor, PrevColor, BlendFactor);
		}

		else 
		{
			o_Color = CurrentColor;
		}
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords);
	}
}

