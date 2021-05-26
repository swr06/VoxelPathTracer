#version 330 core

// Used for reflection denoising
// Based on the assumption that the more rough the surface is, the more denoising is needed

#define INV_SQRT_OF_2PI 0.39894228040143267793994605993439  // 1.0/SQRT_OF_2PI
#define INV_PI 0.31830988618379067153776752674503

out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_PBRTexture;
uniform int u_Radius = 12;

float RoughnessAt = 0.0f;

// https://github.com/BrutPitt/glslSmartDeNoise/ //
vec4 smartDeNoise(sampler2D tex, vec2 uv, float sigma, float kSigma, float threshold)
{
    float radius = round(kSigma*sigma);
    radius = u_Radius;
    float radQ = radius * radius;
    
    float invSigmaQx2 = .5 / (sigma * sigma);      // 1.0 / (sigma^2 * 2.0)
    float invSigmaQx2PI = INV_PI * invSigmaQx2;    // 1.0 / (sqrt(PI) * sigma)
    
    float invThresholdSqx2 = .5 / (threshold * threshold);     // 1.0 / (sigma^2 * 2.0)
    float invThresholdSqrt2PI = INV_SQRT_OF_2PI / threshold;   // 1.0 / (sqrt(2*PI) * sigma)
    
    vec4 centrPx = texture(tex,uv);
    
    float zBuff = 0.0;
    vec4 aBuff = vec4(0.0);
    vec2 size = vec2(textureSize(tex, 0));
    
    for(float x = -radius; x <= radius; x++) 
    {
        float pt = sqrt(radQ - x * x);  

        for(float y = -pt; y <= pt; y++) 
        {
            vec2 d = vec2(x,y);

            float blurFactor = exp(-dot(d , d) * invSigmaQx2) * invSigmaQx2PI; 
            
            vec4 walkPx =  texture(tex,uv+d/size);

            vec4 dC = walkPx - centrPx;
            float deltaFactor = exp(-dot(dC, dC) * invThresholdSqx2) * invThresholdSqrt2PI * blurFactor;
                                 
            zBuff += deltaFactor;
            aBuff += deltaFactor * walkPx;
        }
    }

    return aBuff/zBuff;
}

void main()
{
    RoughnessAt = texture(u_PBRTexture, v_TexCoords).r;
    float EdgeThreshold = mix(0.050f, 0.850f, RoughnessAt * 0.95f);
    o_Color = smartDeNoise(u_Texture, v_TexCoords, 5.0, 2.0, EdgeThreshold).rgb; 
}