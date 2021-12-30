#version 330 core

//#define BE_USELESS

layout (location = 0) out float o_Color;
in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_IntersectionTransversals;

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

float Luminance(vec3 x)
{
	return dot(x, vec3(0.2125f, 0.7154f, 0.0721f));
}

float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

float ShadowSpatial(sampler2D tex, vec2 uv)
{
    const float Diagonal = sqrt(2.0f);

    vec4 CenterPosition = GetPositionAt(u_PositionTexture, v_TexCoords);
    bool Sky = CenterPosition.w < 0.0f;

    vec3 CenterNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).rgb;
    float CenterShadow = texture(u_InputTexture, v_TexCoords).x;
    float BaseLuminance = Luminance(vec3(CenterShadow));
    vec2 TexelSize = 1.0f / vec2(textureSize(tex, 0));

    float Transversal = texture(u_IntersectionTransversals, v_TexCoords).x;
    Transversal = Transversal * 100.0f;

    const float TransversalCutoff = sqrt(2.0f);
    bool SharpShadow = Transversal > 0.0f && Transversal < TransversalCutoff;

    if (SharpShadow||Sky) {
        return CenterShadow;
    }

    bool ReducedKernel = Transversal < TransversalCutoff * 1.414;

    float TotalWeight = 0.0f;
    float TotalShadow = 0.0f;
    
    // 6 * 4 = 24 samples. 

    int XSize = ReducedKernel ? 1 : 3;
    int YSize = ReducedKernel ? 1 : 3;

    float Scale = Transversal > 8.0f ? 1.2f : 1.0f;

    float ClampedTransversal = clamp(Transversal, 0.0f, 10.0f);
    float VarianceEstimate = mix(20.0f, 6.0f, ClampedTransversal / 10.0f)+(Transversal < 6.0f ? 5.0f : 2.0f);

    for (int x = -XSize ; x <= XSize ; x++) {
        for (int y = -YSize ; y <= YSize ; y++) 
        {
            vec2 d = vec2(x, y);
            vec2 SampleCoord = uv + d * TexelSize * 1.1f * Scale;
            vec4 SamplePosition = GetPositionAt(u_PositionTexture, SampleCoord);
            vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).rgb;

            // Sample valiation : 
            if (distance(SamplePosition.xyz, CenterPosition.xyz) > Diagonal ||
                SampleNormal != CenterNormal)
            {
                continue;
            }
            
            float ShadowAt = texture(u_InputTexture, SampleCoord).x;
            float LumaAt = Luminance(vec3(ShadowAt));
            float LuminanceError = clamp(1.0f - clamp(abs(ShadowAt - CenterShadow) / 3.0f, 0.0f, 1.0f), 0.0f, 1.0f);
            float Weight = clamp(pow(LuminanceError, VarianceEstimate), 0.0f, 1.0f); //Kernel * clamp(pow(LuminanceError, 3.5f), 0.0f, 1.0f);
            Weight = clamp(Weight, 0.01f, 1.0f);
            TotalShadow += ShadowAt * Weight;
            TotalWeight += Weight;
        }
    }

    return TotalShadow / max(TotalWeight, 0.01f);;
}

void main()
{
#ifdef BE_USELESS 
    o_Color = texture(u_InputTexture,v_TexCoords).x;  
    return;
#endif


    o_Color = ShadowSpatial(u_InputTexture, v_TexCoords); 
    return;
}