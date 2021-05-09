#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec3 v_RayDirection;
out vec3 v_RayOrigin;
out vec2 v_TexCoords;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform int u_VertCurrentFrame = 0;
uniform vec2 u_VertDimensions = vec2(100000.0f);

void main()
{
	const vec2 bayerSequenceOffsets[16] = vec2[16](vec2(0, 3) / 16.0, vec2(8, 11) / 16.0, vec2(2, 1) / 16.0, vec2(10, 9) / 16.0, vec2(12, 15) / 16.0, vec2(4, 7) / 16.0, vec2(14, 13) / 16.0, vec2(6, 5) / 16.0, vec2(3, 0) / 16.0, vec2(11, 8) / 16.0, vec2(1, 2) / 16.0, vec2(9, 10) / 16.0, vec2(15, 12) / 16.0, vec2(7, 4) / 16.0, vec2(13, 14) / 16.0, vec2(5, 6) / 16.0);

	gl_Position = vec4(a_Position.xy, 1.0f, 1.0f);
	v_TexCoords.xy = a_TexCoords.xy;

	vec2 Position = a_Position;

	vec4 clip = vec4(Position.xy, -1.0, 1.0);

	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	eye.xy += ((bayerSequenceOffsets[int(mod(u_VertCurrentFrame, 12.0f))] * 2.0 - 1.0) / (u_VertDimensions * 2.0f));

	v_RayDirection = vec3(u_InverseView * eye);
	v_RayOrigin = u_InverseView[3].xyz;
}