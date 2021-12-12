#version 330 core

layout (location = 0) out vec3 o_Color;

in vec3 v_RayDirection;

uniform samplerCube u_Cubemap;

void main() {

	o_Color = texture(u_Cubemap, v_RayDirection / length(v_RayDirection)).xyz;

}