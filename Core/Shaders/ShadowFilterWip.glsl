#version 330 core

//#define BE_USELESS

layout (location = 0) out float o_Shadow;
layout (location = 1) out float o_Transversals;

in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_IntersectionTransversals;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform int u_Step;

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

float ShadowHeuristic(vec3 x) {
    return clamp(x.x, 0.000000001f, 1.0f);
}

vec2 ShadowSpatial()
{
	const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

    vec2 TexelSize = 1.0f / vec2(textureSize(u_InputTexture, 0));

    float CenterDepth = texture(u_PositionTexture, v_TexCoords).x;
    vec3 CenterNormal = SampleNormalFromTex(u_NormalTexture, v_TexCoords).rgb;
    float CenterShadow = texture(u_InputTexture, v_TexCoords).x;
    float CenterHeuristic = ShadowHeuristic(vec3(CenterShadow));
    float Transversal = texture(u_IntersectionTransversals, v_TexCoords).x;
    float SampleTransversalCenter = Transversal;
    Transversal = Transversal * 100.0f;

    float LengthPhi = Transversal / 4.0f;
    LengthPhi = LengthPhi * LengthPhi;
    LengthPhi /= 1.5f;

    float TotalWeight = 1.0f;
    float TotalShadow = CenterShadow;
    float TotalTransversals = SampleTransversalCenter;
     
    for (int x = -2 ; x <= 2 ; x++) {
        
        for (int y = -2 ; y <= 2 ; y++) {
            
            if (x == 0 && y == 0) { continue; }

            vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize * float(u_Step);

            float SampleDepth = texture(u_PositionTexture, SampleCoord).x;
            vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
            float SampleShadow = texture(u_InputTexture, SampleCoord).x;
            float SampleHeuristic = ShadowHeuristic(vec3(SampleShadow));
            float SampleTransversal = texture(u_IntersectionTransversals, SampleCoord).x;

            float ShadowWeight = abs(SampleHeuristic - CenterHeuristic) / LengthPhi;
            float DepthWeight = pow(exp(-abs(SampleDepth - CenterDepth)), 10.0f);
            float NormalWeight = pow(max(dot(SampleNormal, CenterNormal), 0.0000000001f), 16.0f);
            float XWeight = AtrousWeights[abs(x)];
			float YWeight = AtrousWeights[abs(y)];

            float TotalCurrentWeight = XWeight * YWeight * NormalWeight * DepthWeight * ShadowWeight;

            TotalShadow += SampleShadow * TotalCurrentWeight;
            TotalTransversals += SampleTransversal * TotalCurrentWeight;
            TotalWeight += TotalCurrentWeight;
        }

    }

    TotalWeight = max(TotalWeight, 0.0000001f);
    TotalShadow /= TotalWeight;
    TotalTransversals /= TotalWeight;

    return vec2(TotalShadow, TotalTransversals);
}

void main()
{
#ifdef BE_USELESS 
    o_Shadow = texture(u_InputTexture,v_TexCoords).x;  
    o_Transversals = texture(u_IntersectionTransversals, v_TexCoords).x;
    return;
#endif


    vec2 Result = ShadowSpatial(); 
    o_Shadow = Result.x;
    o_Transversals = Result.y;
    return;
}