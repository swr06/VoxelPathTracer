#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_PBR;

void main() {
	vec3 Fetch = texture(u_PBR, v_TexCoords).xyz;
	float Smoothness = Fetch.x;
	float Roughness = 1.0f - Smoothness;

	float Metalness = Fetch.y;
	
	float Displacement = 0.0f;

	float AO = 1.0f;

	o_Color = vec4(Roughness, Metalness, Displacement, AO);
}