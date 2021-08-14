#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform int u_LOD = 0;
uniform bool u_HQ = true;

in vec2 v_TexCoords;

bool InThresholdedScreenSpace(vec2 x)
{
	const float b = 0.001f;
    return x.x < 1.0f - b && x.x > b && x.y < 1.0f - b && x.y > b;
}

void main()
{
    vec2 TexelSize = 1.0f / textureSize(u_Texture, 0);
    vec3 TotalBloom = vec3(0.0f); 
    float TotalWeight = 0.0f;

    int KernelSize = 5;

    float GaussianWeights[11] = float[] (0.000003, 0.000229, 0.005977, 0.060598, 0.24173, 0.382925,	0.24173, 0.060598, 0.005977, 0.000229, 0.000003);

    for (int i = -KernelSize; i < KernelSize; i++)
    {
        for (int j = -KernelSize; j < KernelSize; j++)
        {
            vec2 S = v_TexCoords + vec2(i,j) * TexelSize;
            if (!InThresholdedScreenSpace(S)) { continue; }
            float CurrentWeight = GaussianWeights[i + 5] * GaussianWeights[j + 5];
            TotalBloom += textureLod(u_Texture, S, 0.0f).rgb * CurrentWeight;
            TotalWeight += CurrentWeight;
        }
    }

    TotalBloom /= TotalWeight;
    o_Color = pow(TotalBloom, vec3(1.0f / 2.2f));
}