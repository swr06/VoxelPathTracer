#version 330 core

out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InitialTracePositionTexture;


vec3 BilateralUpsample(sampler2D tex, vec2 txc)
{
    vec3 base_normal = texture(u_NormalTexture, txc).rgb;
    float base_depth = texture(u_InitialTracePositionTexture, txc).z;

    const vec2 Kernel[4] = vec2[](
        vec2(0.0f, 1.0f),
        vec2(1.0f, 0.0f),
        vec2(-1.0f, 0.0f),
        vec2(0.0, -1.0f)
    );

    vec2 texel_size = 1.0f / textureSize(tex, 0);

    vec3 color = vec3(0.0f, 0.0f, 0.0f);
    float weight_sum;

    for (int i = 0; i < 4; i++) 
    {
        vec3 sampled_normal = texture(u_NormalTexture, txc + Kernel[i] * texel_size).xyz;
        float nweight = pow(abs(dot(sampled_normal, base_normal)), 32);

        float sampled_depth = texture(u_InitialTracePositionTexture, txc + Kernel[i] * texel_size).z; 
        float dweight = 1.0f / (abs(base_depth - sampled_depth) + 0.001f);

        float computed_weight = nweight * dweight;
        color.rgb += texture(tex, txc + Kernel[i] * texel_size).rgb * computed_weight;
        weight_sum += computed_weight;
    }

    color /= max(weight_sum, 0.2f);
    color = clamp(color, texture(tex, txc).rgb * 0.4f, vec3(1.0f));
    return color;

}

vec3 blur5vertical(sampler2D image, vec2 uv, vec2 resolution) 
{
    vec2 direction = vec2(0.0f, 1.0f);
    vec3 color = vec3(0.0);
    vec2 off1 = vec2(1.3333333333333333) * direction;
    color += texture(image, uv).rgb * 0.29411764705882354;
    color += texture(image, uv + (off1 / resolution)).rgb * 0.35294117647058826;
    color += texture(image, uv - (off1 / resolution)).rgb * 0.35294117647058826;
    return color; 
}
  
vec3 blur5(sampler2D image, vec2 uv) 
{
    vec2 resolution = textureSize(u_Texture, 0);

    vec2 direction = vec2(1.0f, 0.0f);
    vec3 color = vec3(0.0);
    vec2 off1 = vec2(1.3333333333333333) * direction;
    color += blur5vertical(image, uv, resolution) * 0.29411764705882354;
    color += blur5vertical(image, uv + (off1 / resolution), resolution) * 0.35294117647058826;
    color += blur5vertical(image, uv - (off1 / resolution), resolution) * 0.35294117647058826;
    return color; 
}

void main()
{
    o_Color = blur5(u_Texture, v_TexCoords);
}
