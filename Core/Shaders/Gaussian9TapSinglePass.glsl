#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;

in vec2 v_TexCoords;

vec4 blur9vertical(sampler2D image, vec2 uv) 
{
    vec2 resolution = textureSize(u_Texture, 0);
    vec2 direction = vec2(0.0f, 1.0f);

    vec4 color = vec4(0.0);
    vec2 off1 = vec2(1.3846153846) * direction;
    vec2 off2 = vec2(3.2307692308) * direction;
    color += texture(image, uv) * 0.2270270270;
    color += texture(image, uv + (off1 / resolution)) * 0.3162162162;
    color += texture(image, uv - (off1 / resolution)) * 0.3162162162;
    color += texture(image, uv + (off2 / resolution)) * 0.0702702703;
    color += texture(image, uv - (off2 / resolution)) * 0.0702702703;
    return color;
}
  
vec4 blur9(sampler2D image, vec2 uv) 
{
    vec2 resolution = textureSize(u_Texture, 0);
    vec2 direction = vec2(1.0f, 0.0f);

    vec4 color = vec4(0.0);
    vec2 off1 = vec2(1.3846153846) * direction;
    vec2 off2 = vec2(3.2307692308) * direction;
    color += blur9vertical(image, uv) * 0.2270270270;
    color += blur9vertical(image, uv + (off1 / resolution)) * 0.3162162162;
    color += blur9vertical(image, uv - (off1 / resolution)) * 0.3162162162;
    color += blur9vertical(image, uv + (off2 / resolution)) * 0.0702702703;
    color += blur9vertical(image, uv - (off2 / resolution)) * 0.0702702703;
    return color;
}

void main()
{
    o_Color = blur9(u_Texture, v_TexCoords).rgb;
}