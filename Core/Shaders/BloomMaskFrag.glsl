#version 330 core

layout (location = 0) out vec3 o_Color; 

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_EmissiveTexture;

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

void main()
{
    float Emissivity = texture(u_EmissiveTexture, v_TexCoords).w;
    vec3 Fetch = texture(u_Texture, v_TexCoords).rgb;
    float L = GetLuminance(Fetch);

    if (Emissivity > 0.0025f) {
        o_Color = (Emissivity * 32.0f) * Fetch;
    }

   //else if (L > 8.0f) {
   //    o_Color = Fetch;
   //}
}