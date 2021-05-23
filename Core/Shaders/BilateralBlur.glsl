#version 330 core

layout(location = 0) out vec3 out_Color;

in vec2 v_TexCoords;
  
uniform float u_Sharpness;
uniform vec2  u_InvResolutionDirection; // either set x to 1/width or y to 1/height
uniform sampler2D u_ColorTexture;
uniform sampler2D u_PositionTexture;

uniform float KERNEL_RADIUS = 5;

vec3 BlurFunction(vec2 uv, float r, vec3 center_c, float center_d, inout float w_total)
{
    vec3  c = texture(u_ColorTexture, uv).xyz;
    float d = texture(u_PositionTexture, uv).z;
  
    float BlurSigma = float(KERNEL_RADIUS) * 0.5f;
    float BlurFalloff = 1.0f / (2.0f * BlurSigma * BlurSigma);
  
    float ddiff = (d - center_d) * u_Sharpness;
    float w = exp2(-r * r * BlurFalloff - ddiff * ddiff);
    w_total += w;

    return c * w;
}

void main()
{
    vec3 center_c = texture(u_ColorTexture, v_TexCoords).xyz;
    float center_d = texture(u_PositionTexture, v_TexCoords).z;
  
    vec3  c_total = center_c;
    float w_total = 1.0;
  
    for (float r = 1; r <= KERNEL_RADIUS; ++r)
    {
        vec2 uv = v_TexCoords + u_InvResolutionDirection * r;
        c_total += BlurFunction(uv, r, center_c, center_d, w_total);  
    }
  
    for (float r = 1; r <= KERNEL_RADIUS; ++r)
    {
        vec2 uv = v_TexCoords - u_InvResolutionDirection * r;
        c_total += BlurFunction(uv, r, center_c, center_d, w_total);  
    }

    out_Color = c_total/w_total;
}