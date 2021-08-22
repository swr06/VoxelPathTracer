#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec2 v_TexCoords;

uniform int u_VertCurrentFrame = 0;
uniform vec2 u_VertDimensions = vec2(100000.0f);

void main()
{
	gl_Position = vec4(a_Position.xy, 1.0f, 1.0f);
	v_TexCoords.xy = a_TexCoords.xy;
}