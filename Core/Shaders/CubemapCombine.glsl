#version 330 core

layout (location =  0) out vec3 o_Color;

in vec3 v_RayDirection;

uniform samplerCube u_Skybox;
uniform samplerCube u_NebulaLowRes;
uniform float u_SunVisibility;
uniform mat4 u_RotationMatrix;
uniform float u_NebulaStrength;

void main() {
	vec3 Vv = vec3(u_RotationMatrix*vec4(normalize(v_RayDirection),1.));
	vec3 Vvd = normalize(v_RayDirection);
	float Transmittance = (1.0f - u_SunVisibility);
	Transmittance *= Transmittance;
	o_Color = texture(u_Skybox, Vvd).xyz + (texture(u_NebulaLowRes, Vv).xyz * 0.5f * u_NebulaStrength * Transmittance);
}