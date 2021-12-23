#version 330 core

#define PI 3.1415926535

layout (location = 0) out vec2 o_Data;

in vec2 v_TexCoords;

float dot_c(vec3 a, vec3 b)
{
	return max(dot(a, b), 1e-5);
}

float radicalInverse(uint bits) {
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

vec2 hammersley(int i, int N){
    return vec2(float(i)/float(N), radicalInverse(uint(i)));
}

void pixarONB(vec3 n, out vec3 b1, out vec3 b2){
	float sign_ = sign(n.z);
	float a = -1.0 / (sign_ + n.z);
	float b = n.x * n.y * a;
	b1 = vec3(1.0 + sign_ * n.x * n.x * a, sign_ * b, -sign_ * n.x);
	b2 = vec3(b, sign_ + n.y * n.y * a, -n.y);
}

vec3 importanceSampleGGX(vec2 randomHemisphere, vec3 N, float roughness){
    float a = roughness*roughness;
    float phi = 2.0 * PI * randomHemisphere.x;
    float cosTheta = sqrt((1.0 - randomHemisphere.y) / (1.0 + (a*a - 1.0) * randomHemisphere.y));
    float sinTheta = sqrt(1.0 - cosTheta*cosTheta);
    vec3 H = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
    vec3 tangent;
    vec3 bitangent;
    pixarONB(N, tangent, bitangent);
    tangent = normalize(tangent);
    bitangent = normalize(bitangent);
    vec3 sampleDir = tangent * H.x + bitangent * H.y + N * H.z;
    return normalize(sampleDir);
}

float geometry(float cosTheta, float k){
	return (cosTheta)/(cosTheta*(1.0-k)+k);
}

float smithsIBL(float NdotV, float NdotL, float roughness){
    float k = (roughness * roughness) / 2.0;; 
	return geometry(NdotV, k) * geometry(NdotL, k);
}

vec2 IntegrateBRDF(float NdotV, float roughness, int sampleCount)
{
    vec3 V;
    V.x = sqrt(1.0 - NdotV*NdotV);
    V.y = 0.0;
    V.z = NdotV;
    vec2 result = vec2(0);

    vec3 N = vec3(0.0, 0.0, 1.0);
    
    for(int i = 0; i < sampleCount; i++){
    
        vec2 randomHemisphere = hammersley(i, sampleCount);
        vec3 H  = importanceSampleGGX(randomHemisphere, N, roughness);
        vec3 L  = normalize(2.0 * dot(V, H) * H - V);

        float NdotL = max(L.z, 0.0);
        float NdotH = max(H.z, 0.0);
        float VdotH = dot_c(V, H);

        if(NdotL > 0.0){
            float G = smithsIBL(NdotV, NdotL, roughness);
            float G_Vis = (G * VdotH) / (NdotH * NdotV);
            float Fc = pow(1.0 - VdotH, 5.0);

            result.x += (1.0 - Fc) * G_Vis;
            result.y += Fc * G_Vis;
        }
    }

    return result / float(sampleCount);
}


void main() {
    o_Data = IntegrateBRDF(v_TexCoords.x, v_TexCoords.y, 512);
}