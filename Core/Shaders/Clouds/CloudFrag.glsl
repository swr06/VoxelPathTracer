#version 330 core


#define CLOUD_HEIGHT 70 // temp
#define PI 3.14159265359
#define TAU (3.14159265359 * 2.0f)
#define HALF_PI (3.14159265359 * 0.5f)
#define ONE_OVER_PI (1.0f / 3.14159265359f)
#define INVERSE_PI (1.0f / 3.14159265359f)
#define CHECKERBOARDING

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

uniform float u_Coverage;
uniform vec3 u_SunDirection;
uniform float BoxSize;
uniform float u_DetailIntensity;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform vec2 u_Modifiers;

uniform bool u_Checker;
uniform bool u_UseBayer;
uniform vec2 u_WindowDimensions;


uniform bool u_HighQualityClouds;


const vec3 PlayerOrigin = vec3(0,6200,0); 
const float PlanetRadius = 7773; 
const float AtmosphereRadius = 19773; 
const float Size = AtmosphereRadius - PlanetRadius; 

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

vec3 GetSkyColorAt(vec3 rd) 
{
    vec3 unit_direction = normalize(rd);

    float t = 0.5f * (unit_direction.y + 1.0);
    return (1.0 - t) * vec3(1.0, 1.0, 1.0) +  t * vec3(0.5, 0.7, 1.0);
}

float ConvertValueRange(in float v, in vec2 r1, in vec2 r2)
{
	float ret = (((v - r1.x) * (r2.y - r2.x)) / (r1.y - r1.x)) + r2.x;
	return ret;
}

vec2 RayBoxIntersect(vec3 boundsMin, vec3 boundsMax, vec3 rayOrigin, vec3 invRaydir)
{
	vec3 t0 = (boundsMin - rayOrigin) * invRaydir;
	vec3 t1 = (boundsMax - rayOrigin) * invRaydir;
	vec3 tmin = min(t0, t1);
	vec3 tmax = max(t0, t1);
	
	float dstA = max(max(tmin.x, tmin.y), tmin.z);
	float dstB = min(tmax.x, min(tmax.y, tmax.z));
	
	// CASE 1: ray intersects box from outside (0 <= dstA <= dstB)
	// dstA is dst to nearest intersection, dstB dst to far intersection
	
	// CASE 2: ray intersects box from inside (dstA < 0 < dstB) 
	// dstA is the dst to intersection behind the ray, dstB is dst to forward intersection
	
	// CASE 3: ray misses box (dstA > dstB)
	
	float dstToBox = max(0, dstA);
	float dstInsideBox = max(0, dstB - dstToBox);
	return vec2(dstToBox, dstInsideBox);
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

const mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );

vec2 hash( vec2 p ) {
	p = vec2(dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)));
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

// basic fbm
float noise( in vec2 p ) {
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;
	vec2 i = floor(p + (p.x+p.y)*K1);	
    vec2 a = p - i + (i.x+i.y)*K2;
    vec2 o = (a.x>a.y) ? vec2(1.0,0.0) : vec2(0.0,1.0); //vec2 of = 0.5 + 0.5*vec2(sign(a.x-a.y), sign(a.y-a.x));
    vec2 b = a - o + K2;
	vec2 c = a - 1.0 + 2.0*K2;
    vec3 h = max(0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
	vec3 n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
    return dot(n, vec3(70.0));	
}

float fbm(vec2 n) {
	float total = 0.0, amplitude = 0.1;
	for (int i = 0; i < 7; i++) {
		total += noise(n) * amplitude;
		n = m * n;
		amplitude *= 0.4;
	}
	return total;
}

#define DETAIL

float saturate(float x) {
	return clamp(x,0.0f,1.0f);
}

vec3 DecodeCurlNoise(vec3 c)
{
    return (c - 0.5) * 2.0;
}

float SampleDensity(vec3 p)
{
	p.xyz *= 76.0f;
	vec3 windOffset = vec3(u_Time * 12.6942069420f, 0.0f, 0.0f);
    vec3 pos = p + windOffset;
	float SCALE = 1.0f;
    float noiseScale = max(SCALE * 0.0004f, 0.00001f);
    vec4 LowFreqNoise = texture(u_CloudNoise, pos * noiseScale);
    float LowFreqFBM = (LowFreqNoise.g * 0.625f) +
                       (LowFreqNoise.b * 0.25f)  +
                       (LowFreqNoise.a * 0.125f);
	LowFreqFBM = clamp(LowFreqFBM, 0.00000001f, 1.0f);
	float X = 1.0f; // idk
	vec2 ShapeNoiseModifier = vec2(0.8f + u_Modifiers.x, 1.1f + u_Modifiers.y);
	float Sample = remap(LowFreqNoise.r * pow(1.2f - X, 0.1f), LowFreqFBM * ShapeNoiseModifier.x, ShapeNoiseModifier.y, 0.00000001f, 1.0f);
    float cloudCoverage = u_Coverage * 0.7f;
    Sample = saturate(remap(Sample, 0.0f, 1.0f, 0.0f, 1.0f));
    Sample *= cloudCoverage;
	return clamp(Sample*1.1f, 0.0f, 10.0f);
}

float hg(float a, float g) 
{
    float g2 = g * g;
    return (1.0f - g2) / (4 * PI * pow(1.0f + g2 - 2.0f * g * (a), 1.5));
}

vec3 Beer(in vec3 v)
{
	return exp(-v);
}

float Beer (in float v)
{
	return exp(-v);
}

float hgPhase(float x, float g)
{
    float g2 = g * g;
	return 0.25 * ((1.0 - g2) * pow(1.0 + g2 - 2.0*g*x, -1.5));
}

float phase2Lobes(float x)
{
    const float m = 0.6;
    const float gm = 0.8;
    
	float lobe1 = hgPhase(x, 0.8 * gm);
    float lobe2 = hgPhase(x, -0.5 * gm);
    
    return mix(lobe2, lobe1, m);
}

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

int MIN = -2147483648;
int MAX = 2147483647;
int RNG_SEED;

int xorshift(in int value) 
{
    // Xorshift*32
    // Based on George Marsaglia's work: http://www.jstatsoft.org/v08/i14/paper
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    return value;
}

int nextInt(inout int seed) 
{
    seed = xorshift(seed);
    return seed;
}

float nextFloat(inout int seed) 
{
    seed = xorshift(seed);
    // FIXME: This should have been a seed mapped from MIN..MAX to 0..1 instead
    return abs(fract(float(seed) / 3141.592653));
}

float nextFloat(inout int seed, in float max) 
{
    return nextFloat(seed) * max;
}

float nextFloat(inout int seed, in float min, in float max) 
{
    return min + (max - min) * nextFloat(seed);
}

// raymarches up to find the ambient denisty
float RaymarchAmbient(vec3 Point) 
{
	const vec3 Direction = normalize(vec3(0.0f, 1.0f, 0.0f));
	float End = (AtmosphereRadius - Point.y) / Direction.y;  
	End = min(End, 5000.0); 
	vec3 StartPosition = Point; 
	vec3 EndPosition = Point + Direction * End; 
	float PreviousTraversal = 0.0; 
	float ReturnTransmittance = 1.0;
	float Accum = 0.0; 
	float Dither = Bayer16(gl_FragCoord.xy);
	int LightSteps = 8;

	for(int Step = 0; Step < LightSteps; Step++) {

		if(ReturnTransmittance < 0.0001) 
		{
			break; 
		}

		float t = float(Step + Dither) / float(LightSteps); 
		vec3 Position = mix(StartPosition, EndPosition, t); 
		float Traversal = mix(0.0, End, t); 
		float StepSize = Traversal - PreviousTraversal; 
		float Height = Traversal * Direction.y; 
		Accum += SampleDensity(Position * 0.002f) * StepSize; 
		PreviousTraversal = Traversal; 
	}

	return Accum * 10.0f;
}

float PowHalf(int n) {
	return pow(0.5f, float(n));
}

float RaymarchLight(vec3 Point)
{
	vec3 Direction = normalize(u_SunDirection);
	float End = (AtmosphereRadius - Point.y) / Direction.y;  
	End = min(End, 5000.0); 
	vec3 StartPosition = Point; 
	vec3 EndPosition = Point + Direction * End; 
	float x = 0.0; 
	float ReturnTransmittance = 1.0;
	float Accum = 0.0; 
	float Dither = Bayer32(gl_FragCoord.xy);
	int LightSteps = 16;

	for(int Step = 0; Step < LightSteps; Step++) {

		if(ReturnTransmittance < 0.0001) 
			break; 
	
		float t = float(Step + Dither) / float(LightSteps); 
		vec3 Position = mix(StartPosition, EndPosition, t); 
		float Traversal = mix(0.0, End, t); 
		float StepSize = Traversal - x; 

		//Grab the density at this point 
		float Height = Traversal * Direction.y; 
		Accum += SampleDensity(Position * 0.002f) * StepSize; 
		x = Traversal; 
	}

	const float SunAbsorbption = 1.0f;
	//float LightTransmittance = exp(-Accum * SunAbsorbption); 
	return Accum * 16.0f;
}

// Thanks to jess for suggesting this
float ScatterIntegral(float x, float coeff)
{
    float a = -coeff * (1.0 / log(2.0));
    float b = -1.0 / coeff;
    float c =  1.0 / coeff;

    return exp2(a * x) * b + c;
}

float henyey_greenstein_phase_func(float mu)
{
	// Henyey-Greenstein phase function factor [-1, 1]
	// represents the average cosine of the scattered directions
	// 0 is isotropic scattering
	// > 1 is forward scattering, < 1 is backwards
	const float g = 0.76;

	return
	                     (1. - g*g)
	/ //---------------------------------------------
	     ((4. + PI) * pow(1. + g*g - 2.*g*mu, 1.5));
}

float hg2(float a, float g) 
{
      float g2 = g * g;
      return (1.0f - g2) / (4.0f * 3.1415f * pow(1.0f + g2 - 2.0f * g * (a), 1.5f));
}

float powder(float od) 
{
    return 8.0f * (1.0f - 0.97f * exp(-10.0f * od));
}

float henyeyGreensteinPhase(float cosTheta, float g) 
{
	const float norm = 1.0f / TAU;
	float gg = g * g;
	return norm * ((1.0f - gg) / pow(1.0f + gg - 2.0f * g * cosTheta, 3.0f / 2.0f));
}

float CloudPhaseFunction(float cosTheta, in float an, in float od) 
{
    float frontLobe = henyeyGreensteinPhase(cosTheta * an, pow(0.45f, od + 1.0f));
    float backLobe = henyeyGreensteinPhase(cosTheta * an, -pow(0.45f, od + 1.0f));
    float peakLobe = henyeyGreensteinPhase(cosTheta * an, pow(0.9f, od + 1.0f));
    return mix(mix(frontLobe, backLobe, 0.25f), peakLobe, 0.15f);
}

float CloudAmbientPhase(float cosTheta, in float an, in float od) 
{
    float frontLobe = henyeyGreensteinPhase(cosTheta*an, pow(0.35f, od + 1.0f));
    float backLobe  = henyeyGreensteinPhase(cosTheta*an, -pow(0.35f, od + 1.0f));
    return mix(frontLobe, backLobe, 0.5f);
}

// exponential dropoff : 
// 1, 0.5, 0.25, 0.125 etc
const float[8] ScatterKernel = float[8](
				pow(0.5f, float(0)),
				pow(0.5f, float(1)),
				pow(0.5f, float(2)),
				pow(0.5f, float(3)),
				pow(0.5f, float(4)),
				pow(0.5f, float(5)),
				pow(0.5f, float(6)),
				pow(0.5f, float(7))
);

vec3 GetScatter(float DensitySample, float CosTheta, float CosThetaUp, float Phase, vec3 Point, vec3 SunColor, int CurrentStep)
{
	const float SunBrightness = 3.0f;
    float BeersPowderSample = powder(DensitySample);
    float BeersPowder = mix(BeersPowderSample, 1.0, CosTheta * 0.5f + 0.5f);
	float BeersPowderSky = mix(BeersPowderSample, 1.0, CosThetaUp * 0.5f + 0.5f);
	float AmbientMarchResult = RaymarchAmbient(Point);
	float LightMarchResult = RaymarchLight(Point);

	float CloudShadowCoefficient = 0.30585f; // Increase to have a stronger shadow
	vec2 Scatter = vec2(0.0f);

	// fake second scatter :
	for(int ScatterStep = 0; ScatterStep < 8; ScatterStep++)
	{
		float PhaseSun = CloudPhaseFunction(CosTheta, ScatterKernel[ScatterStep], LightMarchResult);
		float PhaseAmbient = CloudAmbientPhase(CosThetaUp, ScatterKernel[ScatterStep], AmbientMarchResult) * PI;
		vec2 S = vec2(PhaseSun * BeersPowder * exp(-LightMarchResult * CloudShadowCoefficient * ScatterKernel[ScatterStep]),
		              PhaseAmbient * BeersPowderSky * exp(-AmbientMarchResult * CloudShadowCoefficient * ScatterKernel[ScatterStep]));
		//vec2 S = vec2(PhaseSun * BeersPowder * exp(-LightMarchResult * CloudShadowCoefficient * ScatterKernel[ScatterStep]),0.0f);
		
		Scatter += S * ScatterKernel[ScatterStep];
   }

	vec3 SunLight = SunColor * Scatter.x;
	float SkyShadow = Scatter.y * (PI * 0.85f);
    return (SunLight*SkyShadow);
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
    vec2 TexelSize = 1.0f / textureSize(u_PositionTex, 0);
    //vec2 TexelSize = 1.0f / u_Dimensions;

    const vec2 Kernel[8] = vec2[8]
	(
		vec2(-1.0, -1.0),
		vec2( 0.0, -1.0),
		vec2( 1.0, -1.0),
		vec2(-1.0,  0.0),
		vec2( 1.0,  0.0),
		vec2(-1.0,  1.0),
		vec2( 0.0,  1.0),
		vec2( 1.0,  1.0)
	);

    for (int i = 0 ; i < 8 ; i++)
    {
        if (texture(u_PositionTex, txc + Kernel[i] * TexelSize).r <= 0.0f)
        {
            return true;
        }
    }

    return false;
}

void main()
{
	o_Data = vec4(0.0f);
	
	// We dont need to ray cast the clouds if there is a hit at the current position
	if (!SampleValid(v_TexCoords))
	{
		return;
	}

	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(fract(u_Time * 120.0f));

	// xorshift once : 
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
	
	vec3 CameraPosition = u_InverseView[3].xyz;
	vec3 Origin = PlayerOrigin + CameraPosition; 
	vec3 Direction = normalize(ComputeRayDirection());
	float Start = (PlanetRadius - Origin.y) / Direction.y;  
	float End = (AtmosphereRadius - Origin.y) / Direction.y;  

	if(Start > End) 
	{
		float temp = End; 
		End = Start; 
		Start = temp; 
	}

	if(End < 0.0)
	{
		return; 
	}

	Start = max(Start, 0.0); 

	vec3 StartPosition = Origin + Direction * Start; 
	vec3 EndPosition = Origin + Direction * End; 
	int StepCount = 16;
	
	float Transmittance = 1.0f;
	float CosAngle = dot(normalize(u_SunDirection), normalize(Direction));
	float Dither;

	if (u_UseBayer)
	{
		int Frame = u_CurrentFrame % 20;
		vec2 BayerIncrement = vec2(Frame * 1.0f, Frame * 0.5f);
		Dither = Bayer256(gl_FragCoord.xy + BayerIncrement);
	}

	else 
	{
		Dither = nextFloat(RNG_SEED);
	}

	//vec3 SkyLight = texture(u_Atmosphere, vec3(g_Direction.x, g_Direction.y, g_Direction.z)).rgb;
	vec3 SkyLight = vec3(0.0f);
	vec3 Scattering = vec3(0.0f);
	vec3 SunColor = vec3(1.0f);
	float x = 0.0f;


	float CosTheta = dot(Direction, normalize(u_SunDirection));
	const vec3 NormalizedUp = normalize(vec3(0.0f, 1.0f, 0.0f));
	float CosThetaUp = dot(Direction, NormalizedUp);

	const float Phase2Lobes = 1.0f;

	for (int i = 0 ; i < StepCount ; i++)
	{
		float t = float(i + Dither) / float(StepCount); 
		float Traversal = mix(0.0, End-Start, t); 
		vec3 CurrentPoint = mix(StartPosition, EndPosition, t); 
		float DensitySample = SampleDensity(CurrentPoint * 0.002f);
		if (DensitySample < 0.001f) {
			continue;
		}

		DensitySample *= 16.0f;
		float StepSize = Traversal - x; 
		Scattering += GetScatter(DensitySample, CosTheta, CosThetaUp, Phase2Lobes, CurrentPoint, SunColor, i) * Transmittance;
		Transmittance *= exp(-DensitySample * StepSize);
		x = Traversal;
	}
	
	Scattering = pow(Scattering, vec3(0.33333333333f));

	// store it!
	o_Data = vec4(Scattering, Transmittance);
}