#version 330 core

#define INV_SQRT_OF_2PI 0.398942280 
#define INV_PI 0.318309886

layout (location = 0) out float o_Color;
in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

vec3 GetRayDirectionAt(vec2 txc)
{
	vec4 clip = vec4(txc * 2.0f - 1.0f, -1.0f, 1.0f);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0f, 0.0f);
	return vec3(u_InverseView * eye);
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(u_InverseView[3].xyz + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

float SmartDenoise(sampler2D tex, vec2 uv, float sigma, float kSigma, float threshold)
{
    float radius;
    radius = 4;
    float radQ = radius * radius;
    float invSigmaQx2 = 0.5f / (sigma * sigma);   
    float invSigmaQx2PI = INV_PI * invSigmaQx2;    
    float invThresholdSqx2 = 0.5f / (threshold * threshold);
    float invThresholdSqrt2PI = INV_SQRT_OF_2PI / threshold;
    float centrPx = texture(tex, uv).r;
    vec4 CenterPosition = GetPositionAt(u_PositionTexture, v_TexCoords);
    vec3 CenterNormal = texture(u_NormalTexture, v_TexCoords).rgb;

    float zBuff = 0.0f;
    float aBuff = 0.0f;
    vec2 size = vec2(textureSize(tex, 0));
    
    for(float x = -radius; x <= radius; x++) 
    {
        float pt = sqrt(radQ - x * x);  

        for(float y = -pt; y <= pt; y++) 
        {
            vec2 d = vec2(x, y);
            vec2 SampleCoord = uv + d / size;
            vec4 SamplePosition = GetPositionAt(u_PositionTexture, SampleCoord);
            vec3 SampleNormal = texture(u_NormalTexture, SampleCoord).rgb;

            // Sample valiation : 
            if (distance(SamplePosition.xyz, CenterPosition.xyz) > 0.6 ||
                SampleNormal != CenterNormal)
            {
                continue;
            }

            // Weights :
            float blurFactor = exp(-dot(d, d) * invSigmaQx2) * invSigmaQx2PI; 
            float walkPx = texture(tex, uv + d / size).r;
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
    o_Color = SmartDenoise(u_InputTexture, v_TexCoords, 5.0, 2.0, 0.1950f); 


    //o_Color = texture(u_InputTexture, v_TexCoords).r;
}