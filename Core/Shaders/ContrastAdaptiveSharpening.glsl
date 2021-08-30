#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform float u_SharpenAmount;

vec3 LinearToSRGB(vec3 x) 
{
    vec3 r;
  
    if(x.x <= 0.0031308f) {
	r.x = 12.92f * x.x;
    } else {
	r.x = 1.055f * pow(x.x, 1.f/2.4f) - 0.055f;
    }

    if(x.y <= 0.0031308f) {
	r.y = 12.92f * x.y;
    } else {
	r.y = 1.055f * pow(x.y, 1.f/2.4f) - 0.055f;
    }

    if(x.z <= 0.0031308f) {
	r.z = 12.92f * x.z;
    } else {
	r.z = 1.055f * pow(x.z, 1.f/2.4f) - 0.055f;
    }
  
    return r;
}

vec3 ContrastAdaptiveSharpening(sampler2D Texture, ivec2 Pixel, float SharpeningAmount)
{
    vec3 a = texelFetch(Texture, Pixel + ivec2(0.0, -1.0), 0).rgb;
    vec3 b = texelFetch(Texture, Pixel + ivec2(-1.0, 0.0), 0).rgb;
    vec3 c = texelFetch(Texture, Pixel + ivec2(0.0, 0.0), 0).rgb;
    vec3 d = texelFetch(Texture, Pixel + ivec2(1.0, 0.0), 0).rgb;
    vec3 e = texelFetch(Texture, Pixel + ivec2(0.0, 1.0), 0).rgb;
    float MinGreen = min(a.g, min(b.g, min(c.g, min(d.g, e.g))));
    float MaxGreen = max(a.g, max(b.g, max(c.g, max(d.g, e.g))));
    float FinalSharpenAmount = sqrt(min(1.0f - MaxGreen, MinGreen) / MaxGreen);
    float w = FinalSharpenAmount * mix(-0.125f, -0.2f, SharpeningAmount);
    return (w * (a + b + d + e) + c) / (4.0f * w + 1.0f);
}

void main() {
    ivec2 Pixel = ivec2(gl_FragCoord.xy);
    vec3 OriginalColor = texelFetch(u_Texture, Pixel, 0).xyz;
    float SharpeningAmount = u_SharpenAmount;
    vec3 SharpenedColor = ContrastAdaptiveSharpening(u_Texture, Pixel, SharpeningAmount);
    o_Color = LinearToSRGB(SharpenedColor);
}
