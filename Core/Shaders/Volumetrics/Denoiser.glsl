// Basic optimized bilateral denoiser with modified gaussian weights

#version 330 core

layout (location = 0) out vec3 o_PointVolumetrics;

in vec2 v_TexCoords;

uniform sampler2D u_InputPointVolumetrics;
uniform sampler2D u_DepthTexture;


vec4 SampleTextureClamped(sampler2D tex, vec2 uv) {
    if (uv.xy == clamp(uv, 0.00001f, 0.99999f)) {
        return texture(tex, uv.xy);
    }

    return vec4(0.0f);
}

float TexelFetchNormalized(vec2 uv) {
    return texelFetch(u_DepthTexture, ivec2(uv * textureSize(u_DepthTexture, 0).xy), 0).x;
}

float GetGBufferWeight(vec2 UV, float base_depth)
{
	return 1.0f;
    return clamp(pow(1.0f / abs(texture(u_DepthTexture,UV).x - base_depth), 1.0f),0.0125f,1.0f);
}

void GaussianVertical(vec2 uv, float BaseDepth, float wg, inout vec3 OutPointVolumetrics) 
{
    vec2 OneOverResolution = 1.0f / textureSize(u_InputPointVolumetrics, 0);
    vec2 OffsetOne = vec2(1.3846153846f) * vec2(0.0f, 1.0f);
    vec2 OffsetTwo = vec2(3.2307692308f) * vec2(0.0f, 1.0f);

    float Weights[5] = float[5](GetGBufferWeight(uv, BaseDepth),                  
                                GetGBufferWeight(uv + (OffsetOne * OneOverResolution), BaseDepth),
                                GetGBufferWeight(uv - (OffsetOne * OneOverResolution), BaseDepth),
                                GetGBufferWeight(uv + (OffsetTwo * OneOverResolution), BaseDepth),
                                GetGBufferWeight(uv - (OffsetTwo * OneOverResolution), BaseDepth));

    vec3 PointVolumetrics = vec3(0.0f);
    vec4 SunVolumetrics = vec4(0.0f);

    PointVolumetrics += SampleTextureClamped(u_InputPointVolumetrics, uv).xyz * 0.2270270270f * Weights[0];
    PointVolumetrics += SampleTextureClamped(u_InputPointVolumetrics, uv + (OffsetOne * OneOverResolution)).xyz * 0.3162162162f  * Weights[0];
    PointVolumetrics += SampleTextureClamped(u_InputPointVolumetrics, uv - (OffsetOne * OneOverResolution)).xyz * 0.3162162162f  * Weights[0];
    PointVolumetrics += SampleTextureClamped(u_InputPointVolumetrics, uv + (OffsetTwo * OneOverResolution)).xyz * 0.0702702703f  * Weights[0];
    PointVolumetrics += SampleTextureClamped(u_InputPointVolumetrics, uv - (OffsetTwo * OneOverResolution)).xyz * 0.0702702703f  * Weights[0];

    // Weight the above weighted sum with the gaussian horizontal weight :
    OutPointVolumetrics += PointVolumetrics * wg; 
}
  
void Denoise(vec2 UV, float BaseDepth, out vec3 PointVolumetrics) 
{
    vec2 OneOverResolution = 1.0f / textureSize(u_InputPointVolumetrics, 0);

    PointVolumetrics = vec3(0.0f);

    vec2 OffsetOne = vec2(1.3846153846f) * vec2(1.0f, 0.0f);
    vec2 OffsetTwo = vec2(3.2307692308f) * vec2(1.0f, 0.0f);

    GaussianVertical(UV, BaseDepth, 0.2270270270f, PointVolumetrics);
    GaussianVertical(UV + (OffsetOne * OneOverResolution), BaseDepth, 0.3162162162f, PointVolumetrics);
    GaussianVertical(UV - (OffsetOne * OneOverResolution), BaseDepth, 0.3162162162f, PointVolumetrics);
    GaussianVertical(UV + (OffsetTwo * OneOverResolution), BaseDepth, 0.0702702703f, PointVolumetrics);
    GaussianVertical(UV - (OffsetTwo * OneOverResolution), BaseDepth, 0.0702702703f, PointVolumetrics);
}

void main() 
{
    float BaseDepth = texture(u_DepthTexture, v_TexCoords).x; // for gbuffer weight
    Denoise(v_TexCoords, BaseDepth, o_PointVolumetrics);
    o_PointVolumetrics = clamp(o_PointVolumetrics, 0.0f, 16.0f);
}