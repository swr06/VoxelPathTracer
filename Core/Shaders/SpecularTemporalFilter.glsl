#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform sampler2D u_PBRTex;

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

bool Valid(vec2 p)
{
	return p.x > 0.0f && p.x < 1.0f && p.y > 0.0f && p.y < 1.0f;
}

vec2 Reprojection(vec3 pos) 
{
	return ProjectPositionPrevious(pos).xy * 0.5f + 0.5f;
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

	minclr -= 0.0025f; 
	maxclr += 0.005f; 
		
	return clamp(texture(u_PreviousColorTexture, reprojected), minclr, maxclr); 
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;

		vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
		CameraOffset *= 0.0f;
		Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);

		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec3 PrevPosition = texture(u_PreviousFramePositionTexture, Reprojected).xyz;
		float d = abs(distance(PrevPosition, CurrentPosition.xyz)); // Disocclusion

		if (Reprojected.x > 0.0 && Reprojected.x < 1.0 && Reprojected.y > 0.0 && Reprojected.y < 1.0 && 
		    d <= 1.5f)
		{
			float RoughnessAt = texture(u_PBRTex, v_TexCoords).r;
			float BlendFactor = 0.8f;
			
			if (RoughnessAt <= 0.2f)
			{
				float RayDistanceAt = CurrentColor.w;

				if (RayDistanceAt > 0.0f)
				{
					float rD = RayDistanceAt * 10.0f;
					vec3 P = CurrentPosition.xyz;
					vec3 V = normalize(u_CurrentCameraPos - P);
					vec3 PotentialEnd = P + V * rD;
			
					vec2 PrevReprojected = Reprojected;
					Reprojected = Reprojection(PotentialEnd);
					
					if (!Valid(Reprojected))
					{
						Reprojected = PrevReprojected;
					}
			
					BlendFactor = 0.99f;
					//o_Color = vec4(1.0f, 0.0f, 0.0f, 1.0f); return;
					//o_Color = vec4(rD / 30.0f); return;
				}
			}

			vec4 PrevColor = GetClampedColor(Reprojected);
			o_Color = mix(CurrentColor, PrevColor, 0.95f);
		}

		else 
		{
			o_Color = CurrentColor.xyzw;
		}
	}

	else 
	{
		o_Color = texture(u_CurrentColorTexture, v_TexCoords);
	}

	o_Color.w = 1.0f;
}

