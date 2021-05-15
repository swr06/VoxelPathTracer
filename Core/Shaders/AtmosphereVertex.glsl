#version 330 core
layout (location = 0) in vec2 a_Position;
layout (location = 0) in vec2 a_TexCoords;

out vec2 v_TexCoords;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InvProjection;
uniform mat4 u_InvView;

out vec3 v_RayDirection;

void main()
{
    v_TexCoords = a_Position;
    gl_Position = vec4(a_Position, 1.0f, 1.0f);

    vec4 clip = vec4(a_Position.xy, -1.0, 1.0);
    vec4 eye = vec4(vec2(u_InvProjection * clip), -1.0, 0.0);
    v_RayDirection = vec3(u_InvView * eye);
}  