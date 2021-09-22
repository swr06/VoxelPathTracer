#version 330 core

//#define USE_PREVIOUS_FRAME

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform int u_CurrentFrame;
uniform sampler2D u_ColorTexture;
uniform sampler2D u_CurrentPositionTexture;
uniform sampler2D u_PreviousColorTexture;

uniform mat4 u_PreviousView;
uniform mat4 u_PreviousProjection;

vec4 SamplePixel(ivec2 px)
{
	return texelFetch(u_ColorTexture, px, 0).rgba;
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
	ivec2 Pixel = ivec2(floor(v_TexCoords*textureSize(u_ColorTexture,0)));

	if (int(gl_FragCoord.x + gl_FragCoord.y) % 2 == CheckerStep)
	{	
		
		const ivec2 Offsets[4] = ivec2[](ivec2(1, 0), ivec2(0, 1), ivec2(-1, 0), ivec2(0, -1));
		

		vec4 Total = vec4(0.0f);
		Total += SamplePixel(Pixel + Offsets[0]);
		Total += SamplePixel(Pixel + Offsets[1]);
		Total += SamplePixel(Pixel + Offsets[2]);
		Total += SamplePixel(Pixel + Offsets[3]);

		Total /= 4.0f;

		o_Color = Total;
	}	

	else 
	{
		o_Color = SamplePixel(Pixel);
	}
}