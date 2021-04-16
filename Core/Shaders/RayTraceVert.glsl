#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec3 v_RayDirection;
out vec3 v_RayOrigin;
out vec2 v_TexCoords;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

void main()
{
	gl_Position = vec4(a_Position.xy, 1.0f, 1.0f);
	v_TexCoords.xy = a_TexCoords.xy;

	vec4 clip = vec4(a_Position.xy, -1.0, 1.0);

	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	v_RayDirection = vec3(u_InverseView * eye);
	v_RayOrigin = u_InverseView[3].xyz;
}