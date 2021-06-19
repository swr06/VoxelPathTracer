#version 330 core

#define CLOUD_HEIGHT 70
#define PI 3.14159265359
#define TAU (3.14159265359 * 2)
#define CHECKERBOARDING
#define DETAIL

#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

layout (location = 0) out vec4 o_Position;
layout (location = 1) out vec3 o_Data;

in vec2 v_TexCoords;

uniform float u_Time;
uniform int u_CurrentFrame;
uniform int u_SliceCount;
uniform vec2 u_Dimensions;

uniform sampler2D u_WorleyNoise;
uniform sampler3D u_CloudNoise;
uniform sampler3D u_CloudDetailedNoise;

uniform sampler2D u_BlueNoise;

uniform float u_Coverage;
uniform vec3 u_SunDirection;
uniform float BoxSize;
uniform float u_DetailIntensity;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform bool u_Checker;
uniform bool u_UseBayer;
uniform vec2 u_WindowDimensions;

const float SunAbsorbption = 0.16f;
const float LightCloudAbsorbption = 2.1f;

vec3 g_Origin;

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
	vec4 sampled_noise;

	vec3 time = vec3(u_Time, 0.0f, u_Time * 0.5f);
	time *= 0.005f;

	sampled_noise = texture(u_CloudNoise, (point.xzy * 0.01f) + time).rgba;

	float perlinWorley = sampled_noise.x;
	vec3 worley = sampled_noise.yzw;
	float wfbm = worley.x * 0.625f + worley.y * 0.125f + worley.z * 0.250f; 
	
	float cloud = remap(perlinWorley, wfbm - 1.0f, 1.0f, 0.0f, 1.0f);
	
	#ifdef DETAIL
	vec4 detail = texture(u_CloudDetailedNoise, (point.xzy * 0.01f) + (time * 2.0f)).rgba;
	float detail_strength = u_DetailIntensity;
	cloud += remap(detail.x, 1.0f - (u_Coverage * 0.6524f), 1.0f, 0.0f, 1.0f) * (0.0354f * detail_strength);
	cloud += remap(detail.y, 1.0f - (u_Coverage * 0.666f), 1.0f, 0.0f, 1.0f) * (0.055f * detail_strength);
	cloud -= remap(detail.z, 1.0f - (u_Coverage * 0.672f), 1.0f, 0.0f, 1.0f) * (0.075f * detail_strength);
	cloud += remap(detail.w, 1.0f - (u_Coverage * 0.684f), 1.0f, 0.0f, 1.0f) * (0.085f * detail_strength);
	#endif

	cloud = remap(cloud, 1.0f - u_Coverage, 1.0f, 0.0f, 1.0f); 

	return clamp(cloud, 0.0f, 2.0f);
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
	int StepCount = 3;
	vec3 ldir = normalize(vec3(u_SunDirection.x, u_SunDirection.y, u_SunDirection.z));

	float tmin, tmax;
	vec3 origin = vec3(g_Origin.x, 0.0f, g_Origin.z);
	vec2 Dist = RayBoxIntersect(origin + vec3(-BoxSize, CLOUD_HEIGHT, -BoxSize), origin + vec3(BoxSize, CLOUD_HEIGHT - 12, BoxSize), p, 1.0f / ldir);
	bool Intersect = !(Dist.y == 0.0f);
	
	if (!Intersect)
	{
		return 1.0f;
	}
	
	tmin = Dist.x;
	tmax = Dist.y;

	float StepSize = tmax / float(StepCount);

	float TotalDensity = 0.0f;
	vec3 CurrentPoint = p + (ldir * StepSize * 0.5f);
	float Dither = nextFloat(RNG_SEED);
	//float Dither = 1.0f;

	for (int i = 0 ; i < StepCount ; i++)
	{
		float DensitySample = SampleDensity(CurrentPoint);
		TotalDensity += max(0.0f, DensitySample * StepSize);
		CurrentPoint += ldir * (StepSize * Dither);
	}

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

float phase2(float a) 
{
	const float x = 0.1;
	const float y = 0.1;
	const float z = 1.2;
	const float w = 0.05;

    float blend = .5;
    float hgBlend = hg(a,x) * (1-blend) + hg(a,-y) * blend;
    return z + hgBlend * w;
}

float RaymarchCloud(vec3 p, vec3 dir, float tmin, float tmax, out float Transmittance, vec3 RayDir)
{
	dir = normalize(dir);
	int StepCount = 10;
	float StepSize = tmax / float(StepCount);

	vec3 CurrentPoint = p + (dir * StepSize * 0.5f);
	float AccumulatedLightEnergy = 0.0f;
	Transmittance = 1.0f;

	float CosAngle = dot(normalize(RayDir), normalize(u_SunDirection));
	float Phase = phase2(CosAngle);
	float Phase2 = henyey_greenstein_phase_func(CosAngle); // todo : check this ? 
	float Phase2Lobes = phase2Lobes(CosAngle);
	Phase = (Phase * 0.7f) + (Phase2 * 1.0f);

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

	for (int i = 0 ; i < StepCount ; i++)
	{
		float DensitySample = SampleDensity(CurrentPoint);
		float BeersPowder = powder(DensitySample);

		//BeersPowder = pow(BeersPowder, 1.75f);
		float Integral = ScatterIntegral(Transmittance, 1.11f);
		float LightMarchSample = RaymarchLight(CurrentPoint);
		AccumulatedLightEnergy += DensitySample * StepSize * LightMarchSample * Transmittance * (Phase * 1.0f);
		Transmittance *= exp(-DensitySample * StepSize * LightCloudAbsorbption);
		CurrentPoint += dir * (StepSize * (Dither));

		if (Transmittance < 0.01f)
		{
			break;
		}
	}
	
	float TotalCloudDensity = AccumulatedLightEnergy;
	return TotalCloudDensity;
}

vec3 ComputeCloudData(in Ray r)
{
	vec3 Output = vec3(0.0f);
	vec3 origin = vec3(g_Origin.x, 0.0f, g_Origin.z);
	vec2 Dist = RayBoxIntersect(origin + vec3(-BoxSize, CLOUD_HEIGHT, -BoxSize), origin + vec3(BoxSize, CLOUD_HEIGHT - 12, BoxSize), r.Origin, 1.0f / r.Direction);
	bool Intersect = !(Dist.y == 0.0f);

	if (Intersect)
	{
		vec3 IntersectionPosition = r.Origin + (r.Direction * Dist.x);
		o_Position.xyz = IntersectionPosition;
		o_Position.w = Dist.y;

		float Transmittance = 1.0f;
		float CloudAt = RaymarchCloud(IntersectionPosition, r.Direction, Dist.x, Dist.y, Transmittance, r.Direction);
		CloudAt = max(CloudAt, 0.0f);
		Transmittance = max(Transmittance, 0.0f);
		Output = vec3(CloudAt, Transmittance, 0.0f);
	}

	return Output;
}

vec3 ComputeRayDirection()
{
	vec2 ScreenSpace = v_TexCoords;
	float TexelSizeX = 1.0f / u_Dimensions.x;
	ScreenSpace.x += (float(int(gl_FragCoord.x + gl_FragCoord.y) % 2 == int(u_CurrentFrame % 2)) * TexelSizeX) * float(u_Checker);
	vec4 Clip = vec4(ScreenSpace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 Eye = vec4(vec2(u_InverseProjection * Clip), -1.0, 0.0);
	vec3 RayDir = vec3(u_InverseView * Eye);

	return RayDir;
}

void main()
{
	o_Position = vec4(0.0f);
	o_Data = vec3(0.0f);

	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(u_Dimensions.x) * int(u_Time * 60.0f);

	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
	
    Ray r;
    r.Origin = u_InverseView[3].xyz;
    r.Direction = normalize(ComputeRayDirection());
	g_Origin = r.Origin;

	vec3 Accumulated = vec3(0.0f);
	const int SAMPLE_COUNT = 1;

	for (int i = 0 ; i < SAMPLE_COUNT; i++)
	{
		Accumulated += ComputeCloudData(r);
	}

	Accumulated = Accumulated * (1.0f / SAMPLE_COUNT);
	o_Data = Accumulated;
}