// Original shader & technique by UglySwedishFish 

#version 330 core

#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform int u_KernelSize;
uniform sampler2D u_InputTexture;
uniform float u_Scale;

void main() {

    vec2 Kernel[4] = vec2[4](vec2(1.5f, 0.0f), vec2(0.0f, 1.5f), vec2(1.0f, 1.0f), vec2(1.0f,-1.0f));
    
    int KernelSize = u_KernelSize;

    vec4 Total = vec4(1e-8);

    vec2 TexelSize = 1.0f / textureSize(u_InputTexture, 0).xy;

    float Hash = bayer32(gl_FragCoord.xy);
    Hash = mix(Hash, 1.0f, 0.75f);
    
    float StepSize = 4.0f * u_Scale;

    for (int Sample = -KernelSize ; Sample <= KernelSize ; Sample++) {

        Sample = Sample == 0 ? 1 : Sample;

        float CurrentStep = float(Sample) * StepSize * Hash;
        float Weight = abs(float(CurrentStep) * (float(CurrentStep)  / 2.0f)) / 10.0f; 
        Weight = 1.0f / Weight;


        for (int Direction = 0 ; Direction < 4 ; Direction++) {

            vec2 SampleCoord = v_TexCoords + (Kernel[Direction] * vec2(CurrentStep)) * TexelSize;

            if (SampleCoord == clamp(SampleCoord, 0.0f, 1.0f)) 
            {
                Total += texture(u_InputTexture, SampleCoord).xyzw * Weight;
            }
        }

    }


    o_Color = Total;
}