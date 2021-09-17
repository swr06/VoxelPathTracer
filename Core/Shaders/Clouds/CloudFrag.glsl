#version 330 core

#define CLOUD_THICKNESS 1750.0
#define CLOUD_HEIGHT 1500.0
#define CLOUD_TOP (CLOUD_HEIGHT + CLOUD_THICKNESS)

#define PI 3.14159265359
#define TAU (3.14159265359 * 2.0f)
#define HALF_PI (3.14159265359 * 0.5f)
#define ONE_OVER_PI (1.0f / 3.14159265359f)
#define INVERSE_PI (1.0f / 3.14159265359f)
#define CHECKERBOARDING
#define LIGHT_LOD 4.0f 
#define BASE_LOD 0.0f

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
uniform vec3 u_DetailParams;

uniform bool u_Checker;
uniform bool u_UseBayer;
uniform vec2 u_WindowDimensions;
uniform float u_TimeScale = 1.0f;


uniform bool u_HighQualityClouds;


vec3 NoiseKernel[6] = vec3[] 
(
	vec3( 0.38051305,  0.92453449, -0.02111345),
	vec3(-0.50625799, -0.03590792, -0.86163418),
	vec3(-0.32509218, -0.94557439,  0.01428793),
	vec3( 0.09026238, -0.27376545,  0.95755165),
	vec3( 0.28128598,  0.42443639, -0.86065785),
	vec3(-0.16852403,  0.14748697,  0.97460106)
);

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
	float dstToBox = max(0, dstA);
	float dstInsideBox = max(0, dstB - dstToBox);
	return vec2(dstToBox, dstInsideBox);
}

float intersetCloudPlane(vec3 r0, vec3 rd) {
    // vechitNormal = -normal*sign(dot(normal, direction));
    vec3 point = vec3(0.0, CLOUD_HEIGHT, 0.0);
    return clamp(dot(point - r0, vec3(0,-1,0)) / dot(rd, vec3(0,-1,0)), -1.0, 9991999.0);
}

float intersetCloudPlaneTop(vec3 r0, vec3 rd) {
    // vechitNormal = -normal*sign(dot(normal, direction));
    vec3 point = vec3(0.0, CLOUD_TOP, 0.0);
    return clamp(dot(point - r0, vec3(0,-1,0)) / dot(rd, vec3(0,-1,0)), -1.0, 9991999.0);
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



float SampleDensity(vec3 p, float lod)
{
	p.xyz *= 345.0f;

	vec3 TimeOffset = vec3(u_Time * u_TimeScale * 15.694206942069420694206942069420f, 0.0f, 0.0f);
    vec3 pos = p + TimeOffset;
	float SCALE = 1.0f;
    float NOISE_SCALE = max(SCALE * 0.0004f, 0.00001f);
    vec4 LowFreqNoise = textureLod(u_CloudNoise, pos * NOISE_SCALE, lod);
    float LowFreqFBM = (LowFreqNoise.g * 0.625f) +
                       (LowFreqNoise.b * 0.25f)  +
                       (LowFreqNoise.a * 0.125f);
	LowFreqFBM = clamp(LowFreqFBM, 0.00000001f, 1.0f);
	float X = 1.0f; // idfk
	vec2 ShapeNoiseModifier = vec2(0.8f + u_Modifiers.x, 1.1f + u_Modifiers.y);
	float Sample = remap(LowFreqNoise.r * pow(1.2f - X, 0.1f), LowFreqFBM * ShapeNoiseModifier.x, ShapeNoiseModifier.y, 0.00000001f, 1.0f);
    float EffectiveCoverage = u_Coverage * 0.7f;
    Sample = saturate(remap(Sample, 0.0f, 1.0f, 0.0f, 1.0f));
    Sample *= EffectiveCoverage;

	// u_DetailParams.y > 0.025 if the detail contribution is enabled
	if (u_DetailParams.y < 0.025f) {
		return Sample;
	}

	// todo : improve this shit
	vec4 DetailNoise = texture(u_CloudDetailedNoise, pos * NOISE_SCALE * 1.66f);
    float HighFreqFBM = (LowFreqNoise.g * 0.625f) +
                       (LowFreqNoise.b * 0.25f)  +
                       (LowFreqNoise.a * 0.125f);
	float OneMinusBaseFBM = 1.0f - LowFreqFBM;
    float ErodeWeight = pow(OneMinusBaseFBM, u_DetailParams.z);
    float FinalDensity = Sample - (1.0f - HighFreqFBM) * ErodeWeight * u_DetailParams.x;
	return clamp(FinalDensity * 1.1f, 0.0f, 10.0f);
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
	vec3 Direction = normalize(vec3(0.0f, 1.0f, 0.0f));
	float Accum = 0.0; 
	float Dither = Bayer8(gl_FragCoord.xy);
	int LightSteps = 4;
	float Increment = 24.0f;
    vec3 RayStep = Direction;
	float StepSize = 1.0f / float(LightSteps) * 12.0f;
    vec3 CurrentPoint = Point + RayStep * 22.0f * Dither;

	for(int Step = 0; Step < LightSteps; Step++)
	{
		Increment *= 1.5f;
		Accum += SampleDensity(CurrentPoint * 0.002f, LIGHT_LOD) * StepSize; 
		CurrentPoint += RayStep * Increment;
	}

	const float SunAbsorbption = 1.0f;
	return Accum * 30.0f;
}

float PowHalf(int n) {
	return pow(0.5f, float(n));
}

float RaymarchLight(vec3 Point)
{
	vec3 Direction = normalize(u_SunDirection);
	float Accum = 0.0; 
	float Dither = Bayer8(gl_FragCoord.xy);
	int LightSteps = 6;
	float Increment = 24.0f;
    vec3 RayStep = Direction;
	float StepSize = 1.0f / float(LightSteps) * 12.0f;
    vec3 CurrentPoint = Point + RayStep * 22.0f * Dither;

	for(int Step = 0; Step < LightSteps; Step++)
	{
		Increment *= 1.5f;
		Accum += SampleDensity(CurrentPoint * 0.002f, LIGHT_LOD) * StepSize; 
		CurrentPoint += RayStep * Increment;
	}

	const float SunAbsorbption = 1.0f;
	//float LightTransmittance = exp(-Accum * SunAbsorbption); 
	return Accum * 42.0f;
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
		float PhaseSun = pow(CloudPhaseFunction(CosTheta, ScatterKernel[ScatterStep], LightMarchResult), 1.0f);
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

const float IsotropicScatter = 0.25f / PI;

void main()
{
	o_Data = vec4(0.0f);
	
	// We dont need to ray cast the clouds if there is a hit at the current position
	if (!SampleValid(v_TexCoords))
	{
		return;
	}

	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(fract(u_Time * 120.0f));

	const int StepCount = u_HighQualityClouds ? 16 : 8;


	// xorshift once : 
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
	
	vec3 CameraPosition = u_InverseView[3].xyz;
	vec3 Direction = normalize(ComputeRayDirection());

	float startDist = intersetCloudPlane(CameraPosition, Direction);
    vec3 startPos = CameraPosition + Direction*startDist;
    float endDist = intersetCloudPlaneTop(CameraPosition, Direction);
    vec3 endPos = CameraPosition + Direction*endDist;
    vec3 rayStep = (endPos - startPos) / float(StepCount);
	int Frame = u_CurrentFrame % 30;
	vec2 BayerIncrement = vec2(Frame * 1.0f, Frame * 0.5f);
	float dither = Bayer32(gl_FragCoord.xy+BayerIncrement);
	vec3 CurrentPoint = startPos + rayStep*dither;
	float StepSize = length(rayStep);


	
	float Transmittance = 1.0f;
	float CosAngle = dot(normalize(u_SunDirection), normalize(Direction));


	//vec3 SkyLight = texture(u_Atmosphere, vec3(g_Direction.x, g_Direction.y, g_Direction.z)).rgb;
	vec3 SkyLight = vec3(0.0f);
	vec3 Scattering = vec3(0.0f);
	vec3 SunColor = vec3(1.0f);


	float CosTheta = dot(Direction, normalize(u_SunDirection));
	const vec3 NormalizedUp = normalize(vec3(0.0f, 1.0f, 0.0f));
	float CosThetaUp = dot(Direction, NormalizedUp);

	const float Phase2Lobes = 1.0f;

	for (int i = 0 ; i < StepCount ; i++)
	{
		float DensitySample = SampleDensity(CurrentPoint * 0.002f, BASE_LOD);
		if (DensitySample < 0.001f) {
			continue;
		}

		DensitySample *= 37.0f;
		Scattering += GetScatter(DensitySample, CosTheta, CosThetaUp, Phase2Lobes, CurrentPoint, SunColor, i) * Transmittance;
		Transmittance *= exp(-DensitySample * StepSize);
		CurrentPoint += rayStep;
	}
	
	Scattering = pow(Scattering, vec3(1.0f / 1.4)); 

	// store it!
	o_Data = vec4(Scattering, Transmittance);
	o_Data.xyz = clamp(o_Data.xyz, 0.0f, 1.0f);
}