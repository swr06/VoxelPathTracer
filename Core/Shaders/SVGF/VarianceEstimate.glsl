#version 330 core

layout (location = 0) out float o_Variance;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_SH;
uniform sampler2D u_CoCg;
uniform sampler2D u_Utility;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

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

// From quake2rtx 
vec3 SHToIrridiance(vec4 shY, vec2 CoCg, vec3 v)
{
    float x = dot(shY.xyz, v);
    float Y = 2.0 * (1.023326f * x + 0.886226f * shY.w);
    Y = max(Y, 0.0);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

bool InScreenSpace(in vec2 v) 
{
    return v.x < 1.0f && v.x > 0.0f && v.y < 1.0f && v.y > 0.0f;
}

void main()
{ 
    vec2 Dimensions = textureSize(u_SH, 0);
    vec2 TexelSize = 1.0f / Dimensions;

    vec3 BasePosition = GetPositionAt(v_TexCoords).xyz;
    vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;

    vec3 BaseUtility = texture(u_Utility, v_TexCoords).xyz;

    vec4 BaseSH = texture(u_SH, v_TexCoords).xyzw;
    vec2 BaseCoCg = texture(u_CoCg, v_TexCoords).xy;
    vec3 BaseRadiance = SHToIrridiance(BaseSH, BaseCoCg, normalize(BaseNormal));
    float BaseLuminosity = GetLuminance(BaseRadiance);

    float SPP = BaseUtility.x;
    float BaseMoment = BaseUtility.y;

    float TotalWeight = 0.0f;
    float TotalMoment = 0.0f;
    vec3 TotalRadiance = vec3(0.0f);

    float Variance = 0.0f;

    if (SPP < 4.0f)
    {
        for (int x = -2 ; x <= 2 ; x++) 
        {
            for (int y = -2 ; y <= 2 ; y++)
            {
                vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;

                if (!InScreenSpace(SampleCoord)) { continue; }

                vec3 SampleNormal = texture(u_NormalTexture, SampleCoord).xyz;
                vec3 SamplePosition = GetPositionAt(SampleCoord).xyz;
                vec3 SampleUtility = texture(u_Utility, SampleCoord).xyz;
                float SampleMoment = SampleUtility.y;

                // Sample SH : 
                vec4 SampleSH = texture(u_SH, SampleCoord);
                vec2 SampleCoCg = texture(u_CoCg, SampleCoord).xy;
                vec3 SampleRadiance = SHToIrridiance(SampleSH, SampleCoCg, SampleNormal);
                float SampleLuminosity = GetLuminance(SampleRadiance);

                // Weights : 
                float PositionError = distance(BasePosition, SamplePosition);
                float PositionWeight = 1.0f / PositionError;
                float NormalWeight = pow(max(dot(BaseNormal, SampleNormal), 0.0f), 12.0f);
                float LuminosityWeight = abs(SampleLuminosity - BaseLuminosity) / 1.0e1;
                float Weight = exp(-LuminosityWeight - PositionWeight - NormalWeight);

                TotalWeight += Weight;
                TotalMoment += SampleMoment * Weight;
                TotalRadiance += SampleRadiance * Weight;
            }
        }

        if (TotalWeight > 0.0f) 
        {
            TotalMoment /= TotalWeight;
            TotalRadiance /= TotalWeight;
        }

        float AccumulatedLuminosity = GetLuminance(TotalRadiance);
        AccumulatedLuminosity = AccumulatedLuminosity * AccumulatedLuminosity;
        Variance = TotalMoment - AccumulatedLuminosity;
	    Variance *= 4.0 / SPP;
    } 


    else 
    {
        Variance = BaseMoment - BaseLuminosity * BaseLuminosity;
    }

    o_Variance = Variance;
}