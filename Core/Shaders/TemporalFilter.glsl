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
uniform bool u_Checkerboard = false; 

uniform int u_CurrentFrame;

bool g_CheckerboardStep;

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

vec3 GetCurrentColor()
{
	vec3 FinalCol;
	
	vec3 CurrentSample = texture(u_CurrentColorTexture, v_TexCoords).rgb;
	
	if (g_CheckerboardStep)
	{
		vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0).xy;
		const vec2 Offsets[4] = vec2[](vec2(1.0f, 0.0f), vec2(0.0f, 1.0f), vec2(-1.0f, 0.0f), vec2(0.0f, -1.0f));
		vec3 Averaged = vec3(0.0f);

		for (int i = 0 ; i < 4 ; i++)
		{
			Averaged += texture(u_CurrentColorTexture, v_TexCoords + (Offsets[i] * TexelSize)).rgb;
		}

		return Averaged / 4.0f;
	}

	return CurrentSample;
}

vec3 GetCurrentColorFast(in vec2 txc)
{
	vec3 CurrentSample = texture(u_CurrentColorTexture, txc).rgb;

	if (g_CheckerboardStep)
	{
		vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0).xy;
		const vec2 Offsets[2] = vec2[](vec2(1.0f, 0.0f), vec2(0.0f, -1.0f));
		vec3 Averaged = vec3(0.0f);

		for (int i = 0 ; i < 2 ; i++)
		{
			Averaged += texture(u_CurrentColorTexture, txc + (Offsets[i] * TexelSize)).rgb;
		}

		return Averaged / 2.0f;
	}

	return CurrentSample;
}

vec3 GetClampedColor(vec2 reprojected)
{
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);

	vec3 minclr = vec3(10000.0f); 
	vec3 maxclr = vec3(-10000.0f); 

	for(int x = -1; x <= 1; x++) 
	{
		for(int y = -1; y <= 1; y++) 
		{
			vec3 Sampled = GetCurrentColorFast(v_TexCoords + (vec2(x,y) * TexelSize)); 
			minclr = min(minclr, Sampled); 
			maxclr = max(maxclr, Sampled); 
		}
	}

	minclr -= 0.075f; 
	maxclr += 0.075f; 
	
	return clamp(texture(u_PreviousColorTexture, reprojected).rgb, minclr, maxclr); 

}

void main()
{
	g_CheckerboardStep = false;

	if (u_Checkerboard)
	{
		int CheckerboardStep = u_CurrentFrame % 2 == 0 ? 1 : 0;
		g_CheckerboardStep = int(gl_FragCoord.x + gl_FragCoord.y) % 2 == CheckerboardStep;
	}

	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;

	TexCoord = v_TexCoords;

	vec2 CurrentCoord = TexCoord;
	vec4 CurrentPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgba;
	vec3 CurrentColor = GetCurrentColor();

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

