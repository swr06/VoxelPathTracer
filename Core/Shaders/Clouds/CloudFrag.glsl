#version 330 core

#define CLOUD_HEIGHT 70
#define PI 3.14159265359
#define TAU (3.14159265359 * 2.0f)
#define HALF_PI (3.14159265359 * 0.5f)
#define ONE_OVER_PI (1.0f / 3.14159265359f)
#define CHECKERBOARDING
//#define DETAIL

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
uniform int u_SliceCount;
uniform vec2 u_Dimensions;

uniform sampler3D u_CloudNoise;
//uniform sampler3D u_CloudDetailedNoise;
uniform samplerCube u_Atmosphere;
uniform sampler2D u_BlueNoise;
uniform sampler2D u_PositionTex;

uniform float u_Coverage;
uniform vec3 u_SunDirection;
uniform float BoxSize;
uniform float u_DetailIntensity;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

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

float SampleDensity(in vec3 point)
{
	point.x += 128.0f;
	point.z += 128.0f;
	vec4 sampled_noise;
	vec3 time = vec3(u_Time, 0.0f, u_Time * 0.5f);
	time *= 0.00400f; 
	sampled_noise = texture(u_CloudNoise, (point.xyz * 0.01759750f) + time).rgba;
	float perlinWorley = sampled_noise.x * 1.0f;
	vec3 worley = sampled_noise.yzw;
	float wfbm = worley.x * 0.625f + worley.y * 0.125f + worley.z * 0.250f; 
	float cloud = remap(perlinWorley, wfbm - 1.0f, 1.0f, 0.0f, 1.0f);
	cloud = remap(cloud, 1.0f - u_Coverage, 1.0f, 0.0f, 1.0f); 
	return clamp(cloud, 0.0f, 100.0f);
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

float RaymarchLight(vec3 p)
{
	int StepCount = u_HighQualityClouds ? 8 : 3;
	vec3 ldir = normalize(vec3(u_SunDirection.x, u_SunDirection.y, u_SunDirection.z));

	float StepSize = 8.0f / float(StepCount);
	float TotalDensity = 0.0f;
	vec3 CurrentPoint = p + (ldir * StepSize * 0.5f);
	float Dither;// = Bayer16(gl_FragCoord.xy);
	Dither = 1.0f;

	for (int i = 0 ; i < StepCount ; i++)
	{
		float DensitySample = SampleDensity(CurrentPoint * 0.002f) * 8.6f;
		TotalDensity += max(0.0f, DensitySample * StepSize);
		CurrentPoint += ldir * (StepSize * Dither);
	}

	const float SunAbsorbption = 1.0f;
	float LightTransmittance = exp(-TotalDensity * SunAbsorbption); 
	return LightTransmittance;
}

float powder(float od)
{
	return 1.0 - exp2(-od * 2.0);
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
      float g2 = g*g;
      return (1-g2) / (4*3.1415*pow(1+g2-2*g*(a), 1.5));
}

// Credits : Robobo1221
vec3 GetScatter(float DensitySample, float Phase, vec3 Point, vec3 SunColor, vec3 SkyLight)
{
	const float SunBrightness = 3.0f;
    float Integral = ScatterIntegral(DensitySample, 1.11);
    float BeersPowder = powder(DensitySample * log(2.0));
	float LightMarchResult = RaymarchLight(Point);
	vec3 SunLight = (SunColor * LightMarchResult * BeersPowder) * Phase * HALF_PI * SunBrightness;
    return (SunLight) * Integral * PI;
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
	int StepCount = u_HighQualityClouds ? 24 : 12;
	
	float Transmittance = 1.0f;
	float CosAngle = dot(normalize(u_SunDirection), normalize(Direction));
	float Phase2Lobes = phase2Lobes(CosAngle) * 0.7f;
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

	for (int i = 0 ; i < StepCount ; i++)
	{
		float t = float(i + Dither) / float(StepCount); 
		float Traversal = mix(0.0, End-Start, t); 
		vec3 CurrentPoint = mix(StartPosition, EndPosition, t); 
		float DensitySample = SampleDensity(CurrentPoint * 0.002f) * 2.0f;
		
		if (DensitySample < 0.001f) {
			continue;
		}

		DensitySample *= 10.75f / 2.0f;
		float StepSize = Traversal - x; 
		Scattering += GetScatter(DensitySample, Phase2Lobes, CurrentPoint, SunColor, SkyLight) * Transmittance;
		Transmittance *= exp(-DensitySample * StepSize);
		x = Traversal;
	}
	
	// This is like the most non-physically based thing on earth but idc xD
	Scattering = pow(Scattering, vec3(1.0f / 16.0f));
	o_Data = vec4(Scattering, Transmittance);
}