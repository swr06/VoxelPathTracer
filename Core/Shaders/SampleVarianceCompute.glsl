#version 330 core

layout (location = 0) out float o_Variance;

in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;

float GetLuminance(vec3 color) 
{
    return dot(color, vec3(0.299, 0.587, 0.114));
}

// Estimates variance based on luminance difference 
float GetVarianceEstimate()
{
    vec2 TexelSize = 1.0f / textureSize(u_InputTexture, 0);
    float AverageLuminance = 0.0f;
    int SamplesTaken = 0;
    vec2 Total = vec2(0.0f);

    for (int i = -2 ; i <= 2 ; i++)
    {
        for (int j = -2 ; j <= 2 ; j++)
        {
            vec2 xy = vec2(i,j);
            vec3 SampleColor = texture(u_InputTexture, v_TexCoords + xy*TexelSize).xyz;
            float l = GetLuminance(SampleColor.xyz);
            Total += vec2(l, l * l);
            SamplesTaken++;
        }
    }

    Total /= max(float(SamplesTaken), 1.0f);
    float TotalVariance = max(0.0, Total.y - Total.x * Total.x);
    return TotalVariance;
}

void main()
{
    o_Variance = GetVarianceEstimate();
}