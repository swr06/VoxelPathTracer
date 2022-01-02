#version 330 core

layout (location = 0) out float o_Variance;

in vec2 v_TexCoords;

uniform sampler2D u_RawShadow;

uniform float u_Scale;

float GetHeuristic(float S) 
{
    return S; 

    // We could also fit S to a curve of our liking 

    // return pow(S, 2.0f);
    // return S * pow(1.-S,4.0);

}

float GetVarianceEstimate()
{
    vec2 TexelSize = 1.0f / textureSize(u_RawShadow, 0);
    float AverageLuminance = 0.0f;
    int SamplesTaken = 0;
    vec2 Total = vec2(0.0f);
    float Scale = u_Scale;

    for (int i = -2 ; i <= 2 ; i++)
    {
        for (int j = -2 ; j <= 2 ; j++)
        {
            vec2 xy = vec2(i,j);
            float Sample = texture(u_RawShadow, v_TexCoords + xy*TexelSize*Scale).x;
            float l = GetHeuristic(Sample);
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