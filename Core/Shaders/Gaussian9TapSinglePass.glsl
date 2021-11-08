#version 330 core

layout (location = 0) out vec4 o_Color;

uniform sampler2D u_Texture;

in vec2 v_TexCoords;

vec4 Sampletexture(sampler2D tex, vec2 uv) {
    return texture(tex, clamp(uv,0.000001f,0.999999f));
}

vec4 blur9vertical(sampler2D image, vec2 uv) 
{
    vec2 resolution = textureSize(u_Texture, 0);
    vec2 direction = vec2(0.0f, 1.0f);

    vec4 color = vec4(0.0);
    vec2 off1 = vec2(1.3846153846) * direction;
    vec2 off2 = vec2(3.2307692308) * direction;
    color += Sampletexture(image, uv) * 0.2270270270;
    color += Sampletexture(image, uv + (off1 / resolution)) * 0.3162162162;
    color += Sampletexture(image, uv - (off1 / resolution)) * 0.3162162162;
    color += Sampletexture(image, uv + (off2 / resolution)) * 0.0702702703;
    color += Sampletexture(image, uv - (off2 / resolution)) * 0.0702702703;
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
    o_Color = blur9(u_Texture, v_TexCoords);
}