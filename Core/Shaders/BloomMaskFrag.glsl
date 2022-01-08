#version 330 core

layout (location = 0) out vec3 o_Color; 

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_EmissiveTexture;

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

bool FloatEqual(float x, float y) {
    return abs(x-y) < 0.05f;
}

void main()
{
    float Emissivity = texture(u_EmissiveTexture, v_TexCoords).x;
    vec3 Fetch = texture(u_Texture, v_TexCoords).rgb;

    if (FloatEqual(Emissivity, -8.0f)) {
        Emissivity = 1.2f;
        Fetch = vec3(3.0f,3.0f,2.0f);
    }

    if (FloatEqual(Emissivity, -16.0f)) {
        Emissivity = 1.4f;
        Fetch = vec3(0.6f,0.6f,1.0f);
    }

    Emissivity = max(Emissivity, 0.00000001f);

    float L = GetLuminance(Fetch);

    if (Emissivity > 0.0025f) {
        o_Color = (Emissivity * 6.4f) * Fetch;
    }

    o_Color = max(o_Color, 0.0000000001f);

   //else if (L > 8.0f) {
   //    o_Color = Fetch;
   //}
}