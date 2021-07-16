#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform int u_LOD = 0;
uniform bool u_HQ = true;

in vec2 v_TexCoords;

vec4 blur5vertical(sampler2D image, vec2 uv, vec2 resolution) 
{
    vec2 direction = vec2(0.0f, 1.0f);
    vec4 color = vec4(0.0);
    vec2 off1 = vec2(1.3333333333333333) * direction;
    color += texture2D(image, uv) * 0.29411764705882354;
    color += texture2D(image, uv + (off1 / resolution)) * 0.35294117647058826;
    color += texture2D(image, uv - (off1 / resolution)) * 0.35294117647058826;
    return color; 
}
  
vec4 blur5(sampler2D image, vec2 uv) 
{
    vec2 resolution = textureSize(u_Texture, 0);

    vec2 direction = vec2(1.0f, 0.0f);
    vec4 color = vec4(0.0);
    vec2 off1 = vec2(1.3333333333333333) * direction;
    color += blur5vertical(image, uv, resolution) * 0.29411764705882354;
    color += blur5vertical(image, uv + (off1 / resolution), resolution) * 0.35294117647058826;
    color += blur5vertical(image, uv - (off1 / resolution), resolution) * 0.35294117647058826;
    return color; 
}

void main()
{
    if (u_HQ) {
        float Scale = 1.5f;
        vec2 TexelSize = 1.0f / textureSize(u_Texture, 0);
        vec3 TotalBloom = vec3(0.0f); 
        float TotalWeight = 0.0f;

         for (int i = -5; i < 5; i++)
         {
            for (int j = -5; j < 5; j++)
            {
                float CurrentWeight = pow(1.0 - length(vec2(i, j)) * 0.125f, 6.0);
                TotalBloom += texture(u_Texture, vec2(i, j) * Scale * TexelSize + v_TexCoords).rgb * CurrentWeight;
                TotalWeight += CurrentWeight;
            }
        }

        TotalBloom /= TotalWeight;
        o_Color = TotalBloom;
    } else {
        o_Color = blur5(u_Texture, v_TexCoords).rgb;
    }
}