#version 330 core
#define EPSILON 0.01f

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

float Luminance(vec3 x) {
    return dot(x, vec3(0.2722287168, 0.6740817658, 0.0536895174));
}

float GetSat(vec3 x) { 
    return length(x); 
}

float CASWeight(vec3 x) {
    //return GetSat(x);
    //return Luminance(x);
    return max(x.g, EPSILON);
}

vec3 ContrastAdaptiveSharpening(sampler2D Texture, ivec2 Pixel, float SharpeningAmount)
{
    // Samples 
    vec3 a = texelFetch(Texture, Pixel + ivec2(0, -1), 0).rgb;
    vec3 b = texelFetch(Texture, Pixel + ivec2(-1, 0), 0).rgb;
    vec3 c = texelFetch(Texture, Pixel + ivec2(0, 0), 0).rgb;
    vec3 d = texelFetch(Texture, Pixel + ivec2(1, 0), 0).rgb;
    vec3 e = texelFetch(Texture, Pixel + ivec2(0, 1), 0).rgb;

    // Weight by luminance 
    float WeightA = CASWeight(a.xyz);
    float WeightB = CASWeight(b.xyz);
    float WeightC = CASWeight(c.xyz);
    float WeightD = CASWeight(d.xyz);
    float WeightE = CASWeight(e.xyz);

    // Calculate bounds :
    float MinWeighter = min(WeightA, min(WeightB, min(WeightC, min(WeightD, WeightE))));
    float MaxWeighter = max(WeightA, max(WeightB, max(WeightC, max(WeightD, WeightE))));

    // Apply weights :
    float FinalSharpenAmount = sqrt(min(1.0f - MaxWeighter, MinWeighter) / MaxWeighter);
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
