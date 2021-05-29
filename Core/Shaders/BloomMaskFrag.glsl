#version 330 core

layout (location = 0) out vec3 o_Color; 

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_EmissiveTexture;

void main()
{
    float Emissivity = texture(u_EmissiveTexture, v_TexCoords).w;
    o_Color = (Emissivity * 12.0f) * texture(u_Texture, v_TexCoords).rgb;
}