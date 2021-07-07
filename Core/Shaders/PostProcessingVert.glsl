#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec2 v_TexCoords;
out vec3 v_RayDirection;
out vec3 v_RayOrigin;

uniform mat4 u_VertInverseView;
uniform mat4 u_VertInverseProjection;

void main()
{
	gl_Position = vec4(a_Position, 0.0f, 1.0f);
	v_TexCoords = a_TexCoords;
	vec2 Position = a_Position;
	vec4 clip = vec4(Position.xy, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_VertInverseProjection * clip), -1.0, 0.0);
	v_RayDirection = vec3(u_VertInverseView * eye);
	v_RayOrigin = u_VertInverseView[3].xyz;
}