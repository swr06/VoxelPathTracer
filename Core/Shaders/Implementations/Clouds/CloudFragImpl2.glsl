#version 330 core

#define CLOUD_THICKNESS 1000.0f
#define CLOUD_HEIGHT 1100.0
#define CLOUD_TOP (CLOUD_HEIGHT + CLOUD_THICKNESS)

#define PI 3.14159265359
#define TAU (3.14159265359 * 2.0f)
#define HALF_PI (3.14159265359 * 0.5f)
#define ONE_OVER_PI (1.0f / 3.14159265359f)
#define INVERSE_PI (1.0f / 3.14159265359f)
#define CHECKERBOARDING

// sample at lower lods for the lighting 
// it does not matter too much.
#define LIGHT_LOD 1.5f 
#define AMBIENT_LOD 1.75f 
#define BASE_LOD 0.0f
//#define DITHER_POS_WITH_CURL_NOISE

#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

layout (location = 0) out vec4 o_Data;

in vec2 v_TexCoords;

uniform float u_Time;
uniform int u_CurrentFrame;
uniform int u_SliceCount; // main noise slice count (256)
uniform vec2 u_Dimensions;

uniform sampler3D u_CloudNoise;
uniform sampler3D u_CloudDetailedNoise;
uniform samplerCube u_Atmosphere;
uniform sampler2D u_BlueNoise;
uniform sampler2D u_PositionTex;
uniform sampler2D u_CurlNoise;
uniform sampler2D u_CloudWeatherMap;
uniform sampler2D u_CirrusClouds;

uniform float u_CirrusStrength;
uniform float u_CirrusScale;

uniform float u_Coverage;
uniform vec3 u_SunDirection;
uniform float BoxSize;
uniform float u_DetailIntensity;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform vec2 u_Modifiers;
uniform vec3 u_DetailParams; // scale, enabled, exponent

uniform bool u_Checker;
uniform bool u_UseBayer;
uniform vec2 u_WindowDimensions;
uniform float u_TimeScale = 1.0f;


uniform bool u_HighQualityClouds;
uniform bool u_CurlNoiseOffset;

const float ACCUM_MULTIPLIER = 50.0f;



// used for dither tests 
const vec3 NoiseKernel[6] = vec3[] 
(
	vec3( 0.38051305,  0.92453449, -0.02111345),
	vec3(-0.50625799, -0.03590792, -0.86163418),
	vec3(-0.32509218, -0.94557439,  0.01428793),
	vec3( 0.09026238, -0.27376545,  0.95755165),
	vec3( 0.28128598,  0.42443639, -0.86065785),
	vec3(-0.16852403,  0.14748697,  0.97460106)
);


// Ray struct 
struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

vec2 g_BayerIncrement; // animate the bayer dither to use the temporal filter more effectively

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

float ConvertValueRange(in float v, in vec2 r1, in vec2 r2)
{
	float ret = (((v - r1.x) * (r2.y - r2.x)) / (r1.y - r1.x)) + r2.x;
	return ret;
}

float RayBasePlaneIntersection(vec3 r0, vec3 rd) 
{
    const vec3 P = vec3(0.0, CLOUD_HEIGHT, 0.0);
    return clamp(dot(P - r0, vec3(0.0f, -1.0f ,0.0f)) / dot(rd, vec3(0.0f, -1.0f, 0.0f)), -1.0f, 9991999.0f);
}

float RayBasePlaneIntersectionTop(vec3 r0, vec3 rd)
{
    const vec3 P = vec3(0.0, CLOUD_TOP, 0.0);
    return clamp(dot(P - r0, vec3(0.0f, -1.0f ,0.0f)) / dot(rd, vec3(0.0f, -1.0f, 0.0f)), -1.0f, 9991999.0f);
}


float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

vec4 SampleNoise(in vec3 p)
{
	vec4 sampled_noise = texture(u_CloudNoise, vec3(p.xzy * 0.01f)).rgba;
	return sampled_noise;
}

#define DETAIL

float saturate(float x) {
	return clamp(x,0.0f,1.0f);
}

vec3 DecodeCurlNoise(vec3 c)
{
    return (c - 0.5) * 2.0;
}

float RemapClamped(float original_value, float original_min, float original_max, float new_min, float new_max)
{
    return new_min + (saturate((original_value - original_min) / (original_max - original_min)) * (new_max - new_min));
}

vec3 SampleWeather(vec3 pos)
{
	return texture(u_CloudWeatherMap, pos.xz).xyz;
}

float phaseg(float x, float g){
    float gg = g * g;
    return (gg * -0.25 /3.14 + 0.25 /3.14) * pow(-2.0 * (g * x) + (gg + 1.0), -1.5);
}

float clamp01(float x) {
	return clamp(x,0.0f,1.0f);
}

float map(float n, float start1, float stop1, float start2, float stop2) {
    return ((n - start1) / (stop1 - start1)) * (stop2 - start2) + start2;
}


float GetHeightFraction(vec3 P) {
	float Altitude = (P.y - CLOUD_HEIGHT) / CLOUD_THICKNESS;
	float Fraction = remap(Altitude, 0.0, 0.4, 0.0, 1.0) * remap(Altitude, 0.6, 1.0, 1.0, 0.0);
	return Fraction;
}

#define STRATUS_GRADIENT vec4(0.0, 0.1, 0.2, 0.3)
#define STRATOCUMULUS_GRADIENT vec4(0.02, 0.2, 0.48, 0.625)
#define CUMULUS_GRADIENT vec4(0.00, 0.1625, 0.88, 0.98)

float getDensityForCloud(float heightFraction, float cloudType)
{
	float stratusFactor = 1.0 - clamp(cloudType * 2.0, 0.0, 1.0);
	float stratoCumulusFactor = 1.0 - abs(cloudType - 0.5) * 2.0;
	float cumulusFactor = clamp(cloudType - 0.5, 0.0, 1.0) * 2.0;
	vec4 baseGradient = stratusFactor * STRATUS_GRADIENT + stratoCumulusFactor * STRATOCUMULUS_GRADIENT + cumulusFactor * CUMULUS_GRADIENT;
	return smoothstep(baseGradient.x, baseGradient.y, heightFraction) - smoothstep(baseGradient.z, baseGradient.w, heightFraction);

}

float SampleDensity(vec3 p, float lod)
{
	vec3 OriginalP=p;
	float OriginalY = p.y;
	p*=0.002f;
	lod = clamp(lod, 0.0f, 1.0f);
	p.xyz *= 425.0f;

	vec3 TimeOffset = vec3((u_Time*u_TimeScale) * 24.694206942069420694206942069420f, 0.0f, (u_Time*u_TimeScale) * 4.694206942069420694206942069420f);
    vec3 pos = p + TimeOffset;
    float NOISE_SCALE = max(1.0f * 0.0004f, 0.00001f);
    vec4 LowFreqNoise = textureLod(u_CloudNoise, pos * NOISE_SCALE, lod);
    float LowFreqFBM = (LowFreqNoise.g * 0.625f) +
                       (LowFreqNoise.b * 0.25f)  +
                       (LowFreqNoise.a * 0.125f);
	
	// trial and error took over from here :')
	float BaseShape = remap(LowFreqNoise.r, (-(1.0f - LowFreqFBM)) * u_Modifiers.x, 1.0f, 0.0f, 1.0f);
	
	float HeightModifier = OriginalP.y / CLOUD_THICKNESS;
	//HeightModifier = 1.0 - HeightModifier;
	HeightModifier = pow(HeightModifier, 3.0f);

	//float density = getDensityForCloud(HeightModifier, 1.0);
	//BaseShape *= (density/HeightModifier);
	
	
	

	vec3 Weather = SampleWeather(p*0.00008f);

	float CloudCoverage = Weather.x / 8.0f;// (pow(HeightModifier, 3.0f)) / 5.0f; // / 16.0f;

	float RemappedToCoverage = remap(BaseShape, CloudCoverage, 1.0, 0.0, 1.0);

	if (u_CurlNoiseOffset) {
		vec3 CurlNoise = DecodeCurlNoise(texture(u_CurlNoise,pos.xz).xyz);
		p += CurlNoise*2.0f;
	}
	
	vec4 DetailNoise = textureLod(u_CloudDetailedNoise, vec3(p * NOISE_SCALE * 9.0f) * u_DetailParams.x, lod).xyzw;
	float HighFreqFBM = (DetailNoise.g * 0.625f) +
                       (DetailNoise.b * 0.25f)  +
                       (DetailNoise.a * 0.125f);

	RemappedToCoverage = RemappedToCoverage - (HighFreqFBM) * pow((1.0f - RemappedToCoverage), u_DetailParams.z);
	RemappedToCoverage = remap(RemappedToCoverage / 0.5f, (HighFreqFBM) * 0.2f, 1.0f, 0.0f, 1.0f);
	return clamp((RemappedToCoverage*u_Coverage), 0.0f, 1.0f);
} 




vec3 NormalizedSUNDIR;

const float LOG2 = log(2.0);
const float InverseLOG2 = 1.0 / LOG2;

float RaymarchLight(vec3 Point)
{
	vec3 Direction = NormalizedSUNDIR;
	float Accum = 0.0; 
	int StepCount = 8;

	float StepSize = (CLOUD_TOP - Point.y) / StepCount;

	for (int i = 0 ; i < StepCount ; i++) {
		float DensityAt = SampleDensity(Point, 0.0f);
		Accum += DensityAt;
		Point += Direction * StepSize;
	}

	float dM = 0.04f;

	return exp2(-Accum * StepSize * dM * InverseLOG2 * 0.125f);
}

float RaymarchAmbient(vec3 Point)
{
	vec3 Direction = vec3(0.0f, 1.0f, 0.0f);
	float Accum = 0.0; 
	int StepCount = 8;

	float StepSize = (CLOUD_TOP - Point.y) / StepCount;

	for (int i = 0 ; i < StepCount ; i++) {
		float DensityAt = SampleDensity(Point, 0.0f);
		Accum += DensityAt;
		Point += Direction * StepSize;
	}

	float dM = 0.04f;

	return exp2(-Accum * StepSize * dM * InverseLOG2 * 0.125f);
}

// Phase functions ->
float phaseg8(float x) 
{
	const float g = 0.8 * 0.7f;
	const float g2 = g * g;
	const float g3 = log2((g2 * -0.25 + 0.25) * 0.8f);
	const float g4 = 1.0 + g2;
	const float g5 = -2.0 * g;
	return exp2(log2(g5 * x + g4) * -1.5 + g3);
}

float phasegm5(float x) 
{
	const float g = -0.5 * 0.7f;
	const float g2 = g * g;
	const float g3 = log2((g2 * -0.25 + 0.25) * (1.0 - 0.8f));
	const float g4 = 1.0 + g2;
    const float g5 = -2.0 * g;
    return exp2(log2(g5 * x + g4) * -1.5 + g3);
}

float phase2lobes(float x) 
{
    return phaseg8(x) + phasegm5(x);
}

float powder(float sampleDensity, float VoL) {
	float powd = 1.0 - exp2(-sampleDensity * 2.0);
	return mix(1.0, powd, clamp01(-VoL * 0.5 + 0.5));
}

vec3 Radiance(vec3 Point, float DensitySample, float Phase, float CosTheta)
{
	vec3 SunColor = vec3(1.0f);
	vec3 SkyColor = vec3(1.0f);

	// 2.0 * PI for multiple scattering 
	// Beer powder term is also accounted for 
	vec3 SunLighting = SunColor * RaymarchLight(Point) * powder(DensitySample, CosTheta) * Phase * (2.0f * PI);

	// Use isotropic phase ->
	vec3 SkyLighting = SkyColor * RaymarchAmbient(Point) * 0.25f * PI; 

	return SunLighting + SkyLighting;
}

vec3 ComputeRayDirection()
{
	// Branchless checkerboarding!
	vec2 ScreenSpace = v_TexCoords;
	float TexelSizeX = 1.0f / u_Dimensions.x;
	ScreenSpace.x += (float(int(gl_FragCoord.x + gl_FragCoord.y) % 2 == int(u_CurrentFrame % 2)) * TexelSizeX) * float(u_Checker);
	vec4 Clip = vec4(ScreenSpace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 Eye = vec4(vec2(u_InverseProjection * Clip), -1.0, 0.0);
	vec3 RayDir = vec3(u_InverseView * Eye);
	return RayDir;
}

bool SampleValid(in vec2 txc)
{
    const ivec2 Kernel[8] = ivec2[8]
	(
		ivec2(-1.0, -1.0),
		ivec2( 0.0, -1.0),
		ivec2( 1.0, -1.0),
		ivec2(-1.0,  0.0),
		ivec2( 1.0,  0.0),
		ivec2(-1.0,  1.0),
		ivec2( 0.0,  1.0),
		ivec2( 1.0,  1.0)
	);

	ivec2 basecoord = ivec2(floor(txc*textureSize(u_PositionTex,0).xy));

    for (int i = 0 ; i < 8 ; i++)
    {
		ivec2 samplecoord = basecoord+Kernel[i];
        //if (texture(u_PositionTex, txc + Kernel[i] * TexelSize * 1.1f).r <= 0.0f)
		if (texelFetch(u_PositionTex, samplecoord,0).x <= 0.0f)
        {
            return true;
        }
    }

    return false;
}

const float IsotropicScatter = 0.25f / PI;

void main()
{
	NormalizedSUNDIR = normalize(u_SunDirection);
	o_Data = vec4(0.0f);
	
	// We dont need to ray cast the clouds if there is a hit at the current position
	if (!SampleValid(v_TexCoords))
	{
		return;
	}

	int Frame = u_CurrentFrame % 30;
	vec2 BayerIncrement = vec2(Frame * 1.0f, Frame * 0.5f); 
	g_BayerIncrement=BayerIncrement;


	const int StepCount = 8;

	vec3 CameraPosition = u_InverseView[3].xyz;
	vec3 Direction = normalize(ComputeRayDirection());

	float T1 = RayBasePlaneIntersection(CameraPosition, Direction);
    vec3 StartPosition = CameraPosition + Direction * T1;
    float T2 = RayBasePlaneIntersectionTop(CameraPosition, Direction);
    vec3 EndPosition = CameraPosition + Direction * T2;
    vec3 RayStep = (EndPosition - StartPosition) / float(StepCount);
	float dither = Bayer64(gl_FragCoord.xy+BayerIncrement);

	vec3 CurrentPoint = StartPosition + RayStep * dither; // * bayer dither
	float StepSize = length(RayStep);
	
	float Transmittance = 1.0f;


	//vec3 SkyLight = texture(u_Atmosphere, vec3(g_Direction.x, g_Direction.y, g_Direction.z)).rgb;
	vec3 Scattering = vec3(0.0f);

	float CosTheta = dot(Direction, NormalizedSUNDIR);
	float PhaseFunction = phase2lobes(CosTheta);


	float rSteps = 1.0f / float(StepCount);

	for (int i = 0 ; i < StepCount ; i++)
	{
		float DensitySample = SampleDensity(CurrentPoint, 0.0) * rSteps * 10.0f;

		if (DensitySample < 0.001f)
		{
			continue;
		}

		//DensitySample *= 1.1f;

		vec3 ScatteringAtPoint = Radiance(CurrentPoint, DensitySample, PhaseFunction, CosTheta);

		float TransmittanceAt = exp2(-DensitySample * InverseLOG2);
		Scattering += ScatteringAtPoint * (-Transmittance * TransmittanceAt + Transmittance);
		Transmittance *= TransmittanceAt;

		CurrentPoint += RayStep;
	}
	

	Scattering = clamp(Scattering, 0.0f, 1.0f);

	// store it!
	o_Data = vec4(Scattering, Transmittance);
	o_Data.xyzw = clamp(o_Data.xyzw, 0.0f, 1.0f);
}