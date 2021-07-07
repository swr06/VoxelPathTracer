#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousFramePositionTexture;

uniform sampler2D u_ReflectionHitData;

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

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

vec4 GetClampedColor(vec2 reprojected, vec3 col)
{
	ivec2 Coord = ivec2(v_TexCoords * Dimensions); 
	vec2 TexelSize = textureSize(u_PreviousColorTexture, 0);

	vec4 minclr = vec4(10000.0f); 
	vec4 maxclr = vec4(-10000.0f); 

	vec2 BestOffset = vec2(0.0f);
	float BestLumaDiff = 10000.0f;

	float BaseLuma = GetLuminance(col);

	for(int x = -2; x <= 2; x++) 
	{
		for(int y = -2; y <= 2; y++) 
		{
			vec4 Fetch = texelFetch(u_CurrentColorTexture, Coord + ivec2(x,y), 0); 
			minclr = min(minclr, Fetch); 
			maxclr = max(maxclr, Fetch); 

			float LumaAt = GetLuminance(Fetch.xyz);
			float LumaDifference = abs(LumaAt - BaseLuma);

			// Find pixel with the nearest luma
			if (LumaDifference < BestLumaDiff)
			{
				BestOffset = vec2(x, y);
				BestLumaDiff = LumaDifference;
			}
		}
	}

	minclr -= 0.0015f; 
	maxclr += 0.0015f; 
		
	return clamp(texture(u_PreviousColorTexture, reprojected + (BestOffset*TexelSize)), minclr, maxclr); 
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(v_RayDirection) * Dist, Dist);
}

void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;

	vec2 CurrentCoord = v_TexCoords;
	vec4 CurrentPosition = GetPositionAt(u_CurrentPositionTexture, v_TexCoords).rgba;

	if (CurrentPosition.a > 0.0f)
	{
		vec2 Reprojected;

		vec3 CameraOffset = u_CurrentCameraPos - u_PrevCameraPos; 
		CameraOffset *= 0.6f;
		Reprojected = Reprojection(CurrentPosition.xyz - CameraOffset);

		float RoughnessAt = texture(u_PBRTex, v_TexCoords).r;
		vec4 CurrentColor = texture(u_CurrentColorTexture, CurrentCoord).rgba;
		vec3 PrevPosition = GetPositionAt(u_PreviousFramePositionTexture, Reprojected).xyz;
		float d = abs(distance(PrevPosition, CurrentPosition.xyz)); // Disocclusion check
		
		if (Reprojected.x > 0.0 && Reprojected.x < 1.0 && Reprojected.y > 0.0 && Reprojected.y < 1.0f && d <= 1.41414f)
		{
			vec4 PrevColor = GetClampedColor(Reprojected, CurrentColor.xyz);
			bool Moved = u_CurrentCameraPos != u_PrevCameraPos;
			float BlendFactor = Moved ? 0.8250f : 0.9250f; 
			o_Color = mix(CurrentColor, PrevColor, BlendFactor);
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

