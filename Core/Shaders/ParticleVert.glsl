#version 330 core

layout (location = 0) in vec3 a_Position;
layout (location = 1) in vec2 a_Txc;
layout (location = 2) in float a_Alpha;
layout (location = 3) in float a_idx;

uniform mat4 u_ViewProjection;
out float v_Alpha;
out float v_Z;
out float v_IDX;
out vec2 v_TexCoords;

void main()
{
	vec3 WorldPosition = a_Position; 
	gl_Position = u_ViewProjection * vec4(WorldPosition, 1.0f);
	v_Alpha = a_Alpha;
	v_Z = gl_Position.z / gl_Position.w;
	v_IDX = a_idx;
	v_TexCoords = a_Txc;
}     