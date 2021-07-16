#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform int u_LOD = 0;
uniform bool u_HQ = true;

in vec2 v_TexCoords;

void main()
{
    float Scale = u_HQ ? 1.5f : 1.7f;
    vec2 TexelSize = 1.0f / textureSize(u_Texture, 0);
    vec3 TotalBloom = vec3(0.0f); 
    float TotalWeight = 0.0f;

    int KernelSize = u_HQ ? 5 : 2;

    for (int i = -KernelSize; i < KernelSize; i++)
    {
        for (int j = -KernelSize; j < KernelSize; j++)
        {
            float CurrentWeight = pow(1.0 - length(vec2(i, j)) * 0.125f, 6.0);
            TotalBloom += texture(u_Texture, vec2(i, j) * Scale * TexelSize + v_TexCoords).rgb * CurrentWeight;
            TotalWeight += CurrentWeight;
        }
    }

    TotalBloom /= TotalWeight;
    o_Color = TotalBloom;
}