#version 330 core

//#define USE_PREVIOUS_FRAME

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform int u_CurrentFrame;
uniform sampler2D u_ColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;

uniform mat4 u_PreviousView;
uniform mat4 u_PreviousProjection;

vec3 SamplePixel(vec2 px)
{
	return texture(u_ColorTexture, px).rgb;
}

vec2 Reprojection(vec3 pos) 
{
	vec3 WorldPos = pos;

	vec4 ProjectedPosition = u_PreviousProjection * u_PreviousView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

void main()
{
	int CheckerStep;
	CheckerStep = u_CurrentFrame % 2;
	vec2 TexelSize = 1.0f / textureSize(u_ColorTexture, 0).xy;

	if (int(gl_FragCoord.x + gl_FragCoord.y) % 2 == CheckerStep)
	{	
		
		const ivec2 Offsets[4] = ivec2[](ivec2(1, 0), ivec2(0, 1), ivec2(-1, 0), ivec2(0, -1));
		
		vec3 Total = vec3(0.0f);
		Total += SamplePixel(v_TexCoords + (Offsets[0] * TexelSize));
		Total += SamplePixel(v_TexCoords + (Offsets[1] * TexelSize));
		Total += SamplePixel(v_TexCoords + (Offsets[2] * TexelSize));
		Total += SamplePixel(v_TexCoords + (Offsets[3] * TexelSize));

		Total /= 4.0f;

		#ifdef USE_PREVIOUS_FRAME
		vec3 FetchedPosition = texture(u_CurrentPositionTexture, v_TexCoords).rgb;
		vec2 Reprojected = Reprojection(FetchedPosition);
		
		if(Reprojected.x > 0.05f && Reprojected.x < 0.95f && Reprojected.y > 0.05f && Reprojected.y < 0.95f)
		{
			vec3 PreviousColor = texture(u_PreviousColorTexture, Reprojected).rgb;
			Total = mix(Total, PreviousColor, 0.6f);
		}
		#endif

		o_Color = Total;
	}	

	else 
	{
		vec3 Fetch = SamplePixel(v_TexCoords);
		o_Color = Fetch;
	}
}