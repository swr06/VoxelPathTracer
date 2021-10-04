#version 330 core


layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;
layout (location = 2) out float o_Variance;


in vec2 v_TexCoords;
in vec3 v_RayOrigin;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_SH;
uniform sampler2D u_CoCg;
uniform sampler2D u_Utility;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform bool DO_SPATIAL = true;
uniform bool AGGRESSIVE_DISOCCLUSION_HANDLING = true;

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 GetPositionAt(vec2 txc)
{
	float Dist = texture(u_PositionTexture, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

float SHToY(vec4 shY)
{
    return max(0, 3.544905f * shY.w);
}

bool InScreenSpace(in vec2 v) 
{
    return v.x < 1.0f && v.x > 0.0f && v.y < 1.0f && v.y > 0.0f;
}

bool CompareFloatNormal(float x, float y) {
    return abs(x - y) < 0.02f;
}

vec3 GetNormalFromID(float n) {
	const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f));
    int idx = int(round(n*10.0f));

    if (idx > 5) {
        return vec3(1.0f, 1.0f, 1.0f);
    }

    return Normals[idx];
}

vec3 SampleNormalFromTex(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}

void main()
{ 
    vec2 Dimensions = textureSize(u_SH, 0);
    vec2 TexelSize = 1.0f / Dimensions;

    vec3 BasePosition = GetPositionAt(v_TexCoords).xyz;
    vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).xyz;

    vec3 BaseUtility = texture(u_Utility, v_TexCoords).xyz;

    vec4 BaseSH = texture(u_SH, v_TexCoords).xyzw;
    vec2 BaseCoCg = texture(u_CoCg, v_TexCoords).xy;
    float BaseLuminosity = SHToY(BaseSH);

    float ACCUMULATED_FRAMES = BaseUtility.x;
    float BaseMoment = BaseUtility.y;

    float TotalWeight = 0.0f;
    float TotalMoment = 0.0f;
    vec4 TotalSH = vec4(0.0f);
    vec2 TotalCoCg = vec2(0.0f);
    float TotalLuminosity = 0.0f;
    float TotalWeight2 = 0.0f;
    float Variance = 0.0f;

    const float FRAME_BIAS = 4.0f;  //u_Framebias////
    const float ACCUMULATED_FRAMES_THRESH = AGGRESSIVE_DISOCCLUSION_HANDLING ? 4.0f+4.0f+FRAME_BIAS : 4.0f+FRAME_BIAS;

    if (ACCUMULATED_FRAMES < ACCUMULATED_FRAMES_THRESH) 
    {
        const float ColorPhi = AGGRESSIVE_DISOCCLUSION_HANDLING ? 5.0f : 5.0f*2.0f; 

        int XKernel = AGGRESSIVE_DISOCCLUSION_HANDLING ? 4 : 1;
        int YKernel = AGGRESSIVE_DISOCCLUSION_HANDLING ? 4 : 1;

        for (int x = -XKernel ; x <= XKernel ; x++) 
        {
            for (int y = -YKernel ; y <= YKernel ; y++)
            {
                vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;

                if (!InScreenSpace(SampleCoord)) { continue; }

                vec3 SamplePosition = GetPositionAt(SampleCoord).xyz;

                // Weights : 
                vec3 PositionDifference = abs(SamplePosition - BasePosition);
                float DistSqr = dot(PositionDifference, PositionDifference);

                if (DistSqr < 1.2f) 
                { 
                    vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
                    vec3 SampleUtility = texture(u_Utility, SampleCoord).xyz;
                    float SampleMoment = SampleUtility.y;

                    // Sample SH : 
                    vec4 SampleSH = texture(u_SH, SampleCoord);
                    vec2 SampleCoCg = texture(u_CoCg, SampleCoord).xy;
                    float SampleLuminosity = SHToY(SampleSH);

                    float NormalWeight = pow(max(dot(BaseNormal, SampleNormal), 0.0f), 16.0f);
                    float LuminosityWeight = abs(SampleLuminosity - BaseLuminosity) / ColorPhi;
                    float Weight = exp(-LuminosityWeight - NormalWeight);
                    float Weight_2 = Weight;

                    Weight = max(Weight, 0.015f);
                    Weight_2 = max(Weight_2, 0.015f);

                    TotalWeight += Weight;
                    TotalMoment += SampleMoment * Weight_2;

                    TotalSH += SampleSH * Weight;
                    TotalCoCg += SampleCoCg * Weight;
                    TotalLuminosity += SampleLuminosity * Weight_2;

                    TotalWeight2 += Weight_2;
                }
            }
        }

        if (TotalWeight > 0.0f) 
        {
            TotalMoment /= TotalWeight2;
            TotalLuminosity /= TotalWeight2;
            TotalCoCg /= TotalWeight;
            TotalSH /= TotalWeight;
        }


        o_SH = TotalSH;
        o_CoCg = TotalCoCg;

        float AccumulatedLuminosity = TotalLuminosity;
        AccumulatedLuminosity = AccumulatedLuminosity * AccumulatedLuminosity;
        Variance = TotalMoment - AccumulatedLuminosity;
	    Variance *= ACCUMULATED_FRAMES_THRESH / ACCUMULATED_FRAMES;
    } 


    else 
    {
        o_SH = BaseSH;
        o_CoCg = BaseCoCg;
        Variance = BaseMoment - BaseLuminosity * BaseLuminosity;
    }

    if (!DO_SPATIAL) {
        o_SH = BaseSH;
        o_CoCg = BaseCoCg;
        Variance = BaseMoment - BaseLuminosity * BaseLuminosity;
    }

    o_Variance = Variance;

    // clamp
    o_SH = clamp(o_SH, -100.0f, 100.0f);
	o_CoCg = clamp(o_CoCg, -10.0f, 100.0f);
	o_Variance = clamp(o_Variance, -1.0f, 50.0f);
}