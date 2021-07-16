#version 330 core

#define INV_SQRT_OF_2PI 0.398942280401432677939946  // 1.0/SQRT_OF_2PI
#define INV_PI 0.31830988618379067153776752674503

layout (location = 0) out float o_Color;
in vec2 v_TexCoords;
uniform sampler2D u_InputTexture;

float SmartDenoise(sampler2D tex, vec2 uv, float sigma, float kSigma, float threshold)
{
    float radius;
    radius = 4;
    float radQ = radius * radius;
    float invSigmaQx2 = 0.5f / (sigma * sigma);      // 1.0 / (sigma^2 * 2.0)
    float invSigmaQx2PI = INV_PI * invSigmaQx2;    // 1.0 / (sqrt(PI) * sigma)
    float invThresholdSqx2 = 0.5f / (threshold * threshold);     // 1.0 / (sigma^2 * 2.0)
    float invThresholdSqrt2PI = INV_SQRT_OF_2PI / threshold;   // 1.0 / (sqrt(2*PI) * sigma)
    float centrPx = texture(tex, uv).r;
    float zBuff = 0.0f;
    float aBuff = 0.0f;
    vec2 size = vec2(textureSize(tex, 0));
    
    for(float x = -radius; x <= radius; x++) 
    {
        float pt = sqrt(radQ - x * x);  

        for(float y = -pt; y <= pt; y++) 
        {
            vec2 d = vec2(x, y);
            float blurFactor = exp(-dot(d, d) * invSigmaQx2) * invSigmaQx2PI; 
            float walkPx =  texture(tex, uv + d / size).r;
            float dC = walkPx - centrPx;
            float deltaFactor = exp(-dot(dC, dC) * invThresholdSqx2) * invThresholdSqrt2PI * blurFactor;
            zBuff += deltaFactor;
            aBuff += deltaFactor * walkPx;
        }
    }

    return aBuff / zBuff;
}

void main()
{
    // High edge threshold because the color delta is HUGE! between areas of shadow and no shadow 
    // Doesnt matter at the end, the hard parts are NOT overblurred
    o_Color = SmartDenoise(u_InputTexture, v_TexCoords, 5.0, 2.0, 0.450000f); 
}