#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec2 v_TexCoords;

uniform int u_VertCurrentFrame;
uniform vec2 u_VertDimensions;

void main()
{
	const vec2 BayerSequenceOffsets[16] = vec2[16](vec2(0, 3) / 16.0, vec2(8, 11) / 16.0, vec2(2, 1) / 16.0, vec2(10, 9) / 16.0, vec2(12, 15) / 16.0, vec2(4, 7) / 16.0, vec2(14, 13) / 16.0, vec2(6, 5) / 16.0, vec2(3, 0) / 16.0, vec2(11, 8) / 16.0, vec2(1, 2) / 16.0, vec2(9, 10) / 16.0, vec2(15, 12) / 16.0, vec2(7, 4) / 16.0, vec2(13, 14) / 16.0, vec2(5, 6) / 16.0);
	gl_Position = vec4(a_Position.xy, 1.0f, 1.0f);
	v_TexCoords.xy = a_TexCoords.xy;
}