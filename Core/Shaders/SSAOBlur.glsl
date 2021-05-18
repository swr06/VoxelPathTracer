#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;

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
    o_Color = blur5(u_Texture, v_TexCoords).rgb;
}