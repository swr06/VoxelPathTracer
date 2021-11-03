#version 330 core
#define EPSILON 0.01f

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform sampler2D u_ColorLUT;
uniform float u_SharpenAmount;
uniform bool u_ColorGrading;
uniform int u_SelectedLUT;
uniform bool u_ColorDither;

float linearToSrgb(float linear){
    float SRGBLo = linear * 12.92;
    float SRGBHi = (pow(abs(linear), 1.0/2.4) * 1.055) - 0.055;
    float SRGB = mix(SRGBHi, SRGBLo, step(linear, 0.0031308));
    return SRGB;
}

float srgbToLinear(float color) {
    float linearRGBLo = color / 12.92;
    float linearRGBHi = pow((color + 0.055) / 1.055, 2.4);
    float linearRGB = mix(linearRGBHi, linearRGBLo, step(color, 0.04045));
    return linearRGB;
}

vec3 linearToSrgb(vec3 linear) {
    vec3 SRGBLo = linear * 12.92;
    vec3 SRGBHi = (pow(abs(linear), vec3(1.0/2.4)) * 1.055) - 0.055;
    vec3 SRGB = mix(SRGBHi, SRGBLo, step(linear, vec3(0.0031308)));
    return SRGB;
}

vec3 srgbToLinear(vec3 color) {
    vec3 linearRGBLo = color / 12.92;
    vec3 linearRGBHi = pow((color + 0.055) / 1.055, vec3(2.4));
    vec3 linearRGB = mix(linearRGBHi, linearRGBLo, step(color, vec3(0.04045)));
    return linearRGB;
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


vec3 Lookup(vec3 color) 
{
    const vec2 InverseSize = vec2(1.0f / 512.0f, 1.0f / 5120.0f);
    mat2 Grid = mat2(vec2(1.0f, InverseSize.y * 512), vec2(0.0f, u_SelectedLUT * InverseSize.y * 512)); 

    // calculate blue component -> used to figure out the quad
    float BlueComponent = color.b * 63.0f;

    // Get Quad 
    vec4 Quad = vec4(0.0f);
    Quad.y = floor(floor(BlueComponent) * 0.125f);
    Quad.x = floor(BlueComponent) - (Quad.y * 8.0f);
    Quad.w = floor(ceil(BlueComponent) * 0.125f);
    Quad.z = ceil(BlueComponent) - (Quad.w * 8.0f);

    // calculate sample pos
    vec4 SamplePosition = (Quad * 0.125f) + (0.123046875f * color.rg).xyxy + 0.0009765625f;

    // fetch
    vec3 Fetch1 = texture2D(u_ColorLUT, SamplePosition.xy * Grid[0] + Grid[1]).rgb;
    vec3 Fetch2 = texture2D(u_ColorLUT, SamplePosition.zw * Grid[0] + Grid[1]).rgb;
    
    // mix
    return mix(Fetch1, Fetch2, fract(BlueComponent));
}

void BasicColorDither(inout vec3 color)
{
    vec3 lestynRGB = vec3(dot(vec2(171.0, 231.0), gl_FragCoord.xy));
    lestynRGB = fract(lestynRGB.rgb / vec3(103.0, 71.0, 97.0));
    color += lestynRGB.rgb / 255.0;
}

void main() {
    ivec2 Pixel = ivec2(gl_FragCoord.xy);
    vec3 OriginalColor = texelFetch(u_Texture, Pixel, 0).xyz;
    float SharpeningAmount = u_SharpenAmount;
    vec3 SharpenedColor = ContrastAdaptiveSharpening(u_Texture, Pixel, SharpeningAmount+0.02f);
    o_Color = linearToSrgb(SharpenedColor);
    o_Color = clamp(o_Color,0.0f,1.0f);

    if (u_ColorGrading) {
        o_Color = Lookup(o_Color);
    }

    if (u_ColorDither) {
		BasicColorDither(o_Color);
    }
}
