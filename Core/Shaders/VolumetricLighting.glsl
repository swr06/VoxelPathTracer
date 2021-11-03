#version 330 core

#define SAMPLE_COUNT 16

#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

layout (location = 0) out float o_Volumetrics;

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_BlueNoiseTexture;

uniform vec3 u_StrongerLightDirection;
uniform vec2 u_Dimensions;

uniform mat4 u_ProjectionMatrix;
uniform mat4 u_ViewMatrix;

vec2 TransformToScreenSpace(vec3 pos)
{
    vec4 ViewSpace = u_ViewMatrix * vec4(pos, 1.0f);
    vec4 Projected = u_ProjectionMatrix * ViewSpace;
    Projected.xyz /= Projected.w;
    Projected.xyz = Projected.xyz * 0.5f + 0.5f;
    return Projected.xy;
} 

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

void main()
{
	vec2 SunScreenSpacePosition = TransformToScreenSpace(u_StrongerLightDirection * 10000.0f); 

	float Density = 0.2f;
    float Weight = 0.1f;
    float Decay = 1.0f;
	float Transmittance = 1.0f;

	float BlueNoiseDither = texture(u_BlueNoiseTexture, v_TexCoords * (u_Dimensions / vec2(256.0f))).r;
	vec2 DeltaTextCoord = vec2(v_TexCoords - SunScreenSpacePosition.xy);
	vec2 StepTexCoord = v_TexCoords.xy;
	DeltaTextCoord *= (1.0 /  float(SAMPLE_COUNT)) * 1.0f;
	DeltaTextCoord *= Bayer32(v_TexCoords * u_Dimensions);

	for(int i = 0; i < SAMPLE_COUNT; i++)
	{
		StepTexCoord -= DeltaTextCoord;

		if (StepTexCoord.x < 0.0f || StepTexCoord.x > 1.0f || StepTexCoord.y < 0.0f || StepTexCoord.y > 1.0f)
		{
			break;
		}

		float samp = float(texture(u_PositionTexture, StepTexCoord).r < 0.0f);
		samp *= Decay * Weight;
		o_Volumetrics += samp;
		Transmittance *= Decay;
	}
}