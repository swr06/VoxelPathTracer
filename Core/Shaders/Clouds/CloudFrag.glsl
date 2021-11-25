#version 330 core

#define PI 3.14159265359
#define pi PI
#define TAU (3.14159265359 * 2.0f)
#define tau TAU
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

in vec2 v_TexCoords;;

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

uniform bool u_LodLighting;

uniform float u_CloudDetailFBMPower;

uniform float u_Coverage;
uniform vec3 u_SunDirection;
uniform vec3 u_ViewerPosition;
uniform float BoxSize;
uniform float u_DetailIntensity;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform vec2 u_Modifiers;
uniform vec3 u_DetailParams; // scale, enabled, exponent

uniform bool u_Checker;
uniform bool u_UseBayer;
uniform vec2 u_WindowDimensions;
uniform vec2 u_JitterValue;
uniform float u_TimeScale = 1.0f;

uniform bool u_HighQualityClouds;
uniform bool u_CurlNoiseOffset;

uniform vec3 u_StepCounts;
uniform bool CHECKER_STEP_COUNT;

uniform float u_SunVisibility;

vec2 g_TexCoords;

// 
const float ACCUM_MULTIPLIER = 50.0f;
const float DENSITY_MULTIPLIER = 50.0f;
//

// Physical values ->

const float EarthRadius = 600000.0f;
const float SphereInnerRadius = 5000.0f;
const float SphereOuterRadius = 17000.0f;
const float PlanetRadius = 6371e3;
const float AtmosphereRadius = 6471e3;

const float CloudAltitudeMin = 1500.0f;
const float CloudThickness = 1600.0f;
const float CloudAltitudeMax = CloudAltitudeMin + CloudThickness;

// Ray struct 
struct Ray
{
	vec3 Origin;
	vec3 Direction;
};


vec2 g_BayerIncrement; // animate the bayer dither to use the temporal filter more effectively

// constants 
const float LOG2 = log(2.0);
const float InverseLOG2 = 1.0 / LOG2;



float clamp01(float x) {
	return clamp(x,0.0f,1.0f);
}

float GetCloudFade(float Distance) {
	return exp2(-Distance * 0.00001f);
}

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


float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

vec4 SampleNoise(in vec3 p)
{
	vec4 sampled_noise = texture(u_CloudNoise, vec3(p.xzy * 0.01f)).rgba;
	return sampled_noise;
}

float MiePhase(float x, float g){
    float gg = g * g;
    float xx = x * x;
    return (3.0 / (8.0 * pi) * (1.0 - gg) / (2.0 + gg)) * (1.0 + xx) / pow(1.0 + gg - 2.0 * g * x, 1.5);
}

float hg_phase(float c, float g)
{
    return ((1.0 - g * g) / pow(1.0 + g * g - 2.0 * g * c, 1.5)) / (4.0 * pi);
}

float hg(float c, float t)
{
    const float g0 = 0.15;
    const float g1 = 0.75;
    float phase0 = hg_phase(c,  pow(g0, t + 1.0));
    float phase1 = hg_phase(c, -pow(g0, t + 1.0));
    return mix(phase0, phase1, g1);
}

float phaseg(float x, float g)
{
    float gg = g * g;
    return (gg * -0.25 /3.14 + 0.25 /3.14) * pow(-2.0 * (g * x) + (gg + 1.0), -1.5);
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



float map(float n, float start1, float stop1, float start2, float stop2) {
    return ((n - start1) / (stop1 - start1)) * (stop2 - start2) + start2;
}

// Found using trial and error 
// Not accurate 
float GetHeightWeight(vec3 P, float Wx) 
{
	float Altitude = (P.y - CloudAltitudeMin) / CloudThickness;
	float Weight = remap(Altitude, 0.0f, 0.4f, 0.0f, 1.0f) * remap(Altitude, 0.6f, 1.0f, 1.0f, 0.0f);
	return saturate(pow(Weight, 1.0f));
}

float GetHeightFraction(vec3 P) 
{
	float Altitude = (P.y - CloudAltitudeMin) / CloudThickness;
	return Altitude; //(length(inPos - sphereCenter) - SPHERE_INNER_RADIUS)/(SPHERE_OUTER_RADIUS - SPHERE_INNER_RADIUS);
}


float SampleDensity(vec3 p, float lod)
{
	lod = clamp(lod, 0.0f, 5.0f);

	vec3 OriginalP = p;
	float HeightFraction = GetHeightFraction(p);


	p.xyz += vec3(100.0f, 10.0f, 250.0f);
	p.xyz *= 0.001f;
	p.xyz *= 425.0f;

	//float CloudDensity = getDensityForCloud(HeightFraction, 1.0f);
	//float DensityHeightRatio = (CloudDensity / HeightFraction);

	float Ts = 0.8f;

	vec3 TimeOffset = vec3((u_Time*u_TimeScale*Ts) * 30.69f, 0.0f, (u_Time*u_TimeScale*Ts) * 4.420f);
    vec3 pos = p + TimeOffset;

    float NOISE_SCALE = max(1.0f * 0.0004f, 0.00001f);
    
	vec4 LowFreqNoise = textureLod(u_CloudNoise, pos * NOISE_SCALE * 0.64f, lod);

	// Construct base shape fbm ->
    float LowFreqFBM = (LowFreqNoise.g * 0.625f) +
                       (LowFreqNoise.b * 0.25f)  +
                       (LowFreqNoise.a * 0.125f);

	LowFreqFBM = pow(LowFreqFBM, 1.0f);
	

	float BaseShape = remap(LowFreqNoise.r, (-(1.0f - LowFreqFBM)) * u_Modifiers.x, 1.0f, 0.0f, 1.0f);
	

	// Sample weather map ->

	vec3 Weather = SampleWeather(p * 0.005f);

	// Erode base shape with height gradient ->

	float HeightWeight = GetHeightWeight(OriginalP, Weather.x);
	BaseShape *= HeightWeight ;


	float CloudCoverage = Weather.x / 4.0f;

	float RemappedToCoverage = remap(BaseShape, CloudCoverage, 1.0, 0.0, 1.0);

	if (u_CurlNoiseOffset) {
		vec3 CurlNoise = DecodeCurlNoise(texture(u_CurlNoise,pos.xz).xyz);

		// Offset by decoded curl noise ->
		p += CurlNoise * 5.0f;
	}
	
	vec4 DetailNoise = textureLod(u_CloudDetailedNoise, vec3(p * NOISE_SCALE * 9.0f) * u_DetailParams.x, clamp(lod-1,0.0f,8.0)).xyzw;

	// construct high freq fbm ->
	float HighFreqFBM = (DetailNoise.g * 0.625f) +
                       (DetailNoise.b * 0.25f)  +
                       (DetailNoise.a * 0.125f);

	HighFreqFBM = pow(HighFreqFBM, 1.2f * u_CloudDetailFBMPower);

	RemappedToCoverage = RemappedToCoverage - (HighFreqFBM) * pow((1.0f - RemappedToCoverage), u_DetailParams.z);
	RemappedToCoverage = remap(RemappedToCoverage / 0.5f, (HighFreqFBM) * 0.2f, 1.0f, 0.0f, 1.0f);
	
	// Multiply by coverage ->
	return clamp((RemappedToCoverage * u_Coverage), 0.0f, 1.0f);
} 


// Normalized sun direction ->
vec3 NormalizedSUNDIR;


float RaymarchLight(vec3 Point)
{
	vec3 Direction = NormalizedSUNDIR;
	float RayLength = 22.0f;
    float Accum = 0.0f;
	int StepCount = clamp(int(u_StepCounts.y), 2, 32);

	Point += Bayer8(gl_FragCoord.xy)*Direction;

    for(int i = 0; i < StepCount; i++)
	{
        Accum += RayLength * SampleDensity(Point, u_LodLighting ? 3.0f : 1.0f);

		// Exponential step size ->
        RayLength *= 1.5; 

		Point += Direction * RayLength;
    }

    return Accum;
}

float RaymarchAmbient(vec3 Point)
{
	const vec3 Direction = normalize(vec3(0.0f, 1.0f, 0.0f));
	float RayLength = 22.0f;
    float Accum = 0.0f;
	int StepCount = clamp(int(u_StepCounts.z),2, 32);

	Point += Bayer4(gl_FragCoord.xy)*Direction;

    for(int i = 0; i < StepCount; i++)
	{
        Accum += RayLength * SampleDensity(Point, u_LodLighting ? 3.0f : 1.0f);

		// Exponential step size ->
        RayLength *= 1.5;

		Point += Direction * RayLength;
    }

    return Accum;
}

float RaymarchBounced(vec3 Point)
{
	const vec3 Direction = normalize(vec3(0.0f, -1.0f, 0.0f));
	float RayLength = 22.0f;
    float Accum = 0.0f;
	int StepCount = 2;

	Point += Bayer4(gl_FragCoord.xy)*Direction;

    for(int i = 0; i < StepCount; i++)
	{
        Accum += RayLength * SampleDensity(Point, u_LodLighting ? 3.0f : 1.0);

		// Exponential step size ->
        RayLength *= 1.5;

		Point += Direction * RayLength;
    }

    return Accum;
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

float powder2(float sampleDensity, float VoL) {
	float powd = 1.0 - exp2(-sampleDensity * 2.0);
	return mix(1.0, powd, clamp01(-VoL * 0.5 + 0.5));
}

vec3 RadianceBasic(vec3 Point, float DensitySample, float Phase, float CosTheta)
{
	vec3 SunColor = vec3(1.0f);
	vec3 SkyColor = vec3(1.0f);

	// 2.0 * PI for multiple scattering 
	// Beer powder term is also accounted for ->
	vec3 SunLighting = SunColor * RaymarchLight(Point) * powder2(DensitySample, CosTheta) * Phase * (2.0f * PI);

	// Use isotropic phase ->
	vec3 SkyLighting = SkyColor * RaymarchAmbient(Point) * 0.25f * PI; 

	return SunLighting + SkyLighting;
}

// Powder terms ->
float Powder(float d) {
	return 1.0f - 0.97f * exp(-10.0f * d);
}

float PowderSecondary(float d) {
	return 1.0f - exp(-2.0f * d);
}

// Powder terms ->
float PowderDirect(float d, float costheta) 
{
	const float G = 0.75f;
	float MiePhase = MiePhase(costheta, G);
	float Powder = PowderSecondary(d);
	Powder = mix(Powder, 1.5f, MiePhase);
	Powder = mix(Powder, 1.0f, costheta*0.5f+0.5f)*Powder;
	return Powder;
}

float PowderIndirect(float d, float costheta) 
{
	float Powder = PowderSecondary(d);
	Powder = mix(Powder, 1.0f, costheta*0.5f+0.5f)*Powder;
	return Powder;
}

// Scatter kernel 
const vec3 ScatterKernel[8] = vec3[8](
	  vec3(pow(0.6f, 0.0f), pow(0.3f, 0.0f), pow(0.8f, 0.0f)),
	  vec3(pow(0.6f, 1.0f), pow(0.3f, 1.0f), pow(0.8f, 1.0f)),
	  vec3(pow(0.6f, 2.0f), pow(0.3f, 2.0f), pow(0.8f, 2.0f)),
	  vec3(pow(0.6f, 3.0f), pow(0.3f, 3.0f), pow(0.8f, 3.0f)),
	  vec3(pow(0.6f, 4.0f), pow(0.3f, 4.0f), pow(0.8f, 4.0f)),
	  vec3(pow(0.6f, 5.0f), pow(0.3f, 5.0f), pow(0.8f, 5.0f)),
	  vec3(pow(0.6f, 6.0f), pow(0.3f, 6.0f), pow(0.8f, 6.0f)),
	  vec3(pow(0.6f, 7.0f), pow(0.3f, 7.0f), pow(0.8f, 7.0f))
);


float PhaseMultiplier = 1.0f;

vec3 ComputeScattering(vec3 Point, float DensitySample, float SampleTransmittance, float TotalTransmittance, vec3 Costheta, vec3 MarchTransmittance)
{
	vec3 PowderPhaseTerms = vec3(PowderDirect(MarchTransmittance.x, Costheta.x), PowderIndirect(MarchTransmittance.y, Costheta.y), PowderIndirect(MarchTransmittance.z, Costheta.z));

	float OneMinusTransmittance = 1.0f - SampleTransmittance;
	float BasePowder = Powder(DensitySample);

	vec3 IntegratedScattering = vec3(0.0f);

	// Account for multiple scattering ->
	// Artificially lower the extinction coefficient (Used a summation of multiple scales)

	for (int ScatterStep = 0 ; ScatterStep < 8 ; ScatterStep++) 
	{
		// x -> Attenuation 
		// y -> Contribution / Absorption
		// z -> Eccentricity extinction 
		
		vec3 Kernel = ScatterKernel[ScatterStep]; 

		// Compute transmittance at each scatter step ->
		float sTransmittance = TotalTransmittance * (OneMinusTransmittance * Kernel.x);
		
		// Weight transmittance by scatter weights ->
		float dTransmittance = sTransmittance * exp(-MarchTransmittance.x * Kernel.y);
		float idTransmittance = sTransmittance * exp(-MarchTransmittance.y * Kernel.y);
		float sidTransmittance = sTransmittance * exp(-MarchTransmittance.z);

		// henyey greenstein phase ->
		vec2 PhaseFunctions;
		PhaseFunctions.x = hg(Costheta.x * Kernel.z, MarchTransmittance.x) * PhaseMultiplier;
		PhaseFunctions.y = hg(Costheta.y * Kernel.z, MarchTransmittance.y) * 2.4f * max(1.0f, PhaseMultiplier * 0.5f);

		// Integrate scattering -> 
		// Transmittance * Powder * Phase (henyey greenstein) 

		IntegratedScattering.x += dTransmittance * PowderPhaseTerms.x * PhaseFunctions.x;
		IntegratedScattering.y += idTransmittance * PowderPhaseTerms.y * PhaseFunctions.y;
		IntegratedScattering.z += sidTransmittance * BasePowder;
	}

	return IntegratedScattering;
}

vec3 ComputeRayDirection()
{
	vec2 ScreenSpace = g_TexCoords;
	//ScreenSpace += (u_JitterValue*4.0f) / u_Dimensions;
	float TexelSizeX = 1.0f / u_Dimensions.x;
	ScreenSpace.x += (float(int(gl_FragCoord.x + gl_FragCoord.y) % 2 == int(u_CurrentFrame % 2)) * TexelSizeX) * float(u_Checker);
	vec4 Clip = vec4(ScreenSpace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 Eye = vec4(vec2(u_InverseProjection * Clip), -1.0, 0.0);
	vec3 RayDir = vec3(u_InverseView * Eye);
	return RayDir;
}

bool SampleValid(in vec2 txc)
{
	ivec2 basecoord = ivec2(floor(txc*textureSize(u_PositionTex,0).xy));

    for (int x = -1 ; x <= 1 ; x++)
    {
		for (int y = -2 ; y <= 2 ; y++) {
		
			ivec2 Kernel = ivec2(x,y);

			ivec2 samplecoord = basecoord+Kernel;
			if (texelFetch(u_PositionTex, samplecoord,0).x <= 0.0f)
			{
				return true;
			}
		}
    }

    return false;
}

vec2 RSI(vec3 origin, vec3 dir, float radius)
{
	float B = dot(origin, dir);
	float C = dot(origin, origin) - radius * radius;
	float D = B * B - C;

	vec2 intersection;

	if (D < 0.0)
	{
		intersection = vec2(-1.0, -1.0);
	} 
	
	else
	{
		D = sqrt(D);
		intersection = -B + vec2(-D, D); 
	}

	return intersection;
}

// Cirrus clouds ->


// By iq ->
// Gives much better quality at effectively no extra cost.
vec4 BetterTexture(sampler2D samp, vec2 uv) 
{
    vec2 textureResolution = textureSize(samp, 0).xy;
	uv = uv*textureResolution + 0.5;
	vec2 iuv = floor(uv);
	vec2 fuv = fract(uv);
	uv = iuv + fuv*fuv*(3.0-2.0*fuv); 
	uv = (uv - 0.5)/textureResolution;
	return texture(samp, uv).xyzw;
}


float GetCirrusDensityAt(vec3 P)
{
	vec2 SamplePosition = P.xz * (0.001f / 2.0f) * u_CirrusScale;
	const float ts = 0.0336f;

	SamplePosition += vec2(u_Time * 0.16666f * ts, u_Time * 0.5f * ts).yx;
	SamplePosition /= 1.4f;

	// Curl Noise offset ->

	vec3 CurlNoise = DecodeCurlNoise(BetterTexture(u_CurlNoise,vec2(SamplePosition.x*2.2f,SamplePosition.y*1.7f)).xyz);
	
	// Offset ->
	SamplePosition += clamp(pow(CurlNoise.xy, vec2(4.0f)) * 0.1f * 0.5 * 1.1f, -8.0f, 8.0f);


	// Shape noise ->
	float CirrusDensity = BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*1.4f,SamplePosition.y*0.92525f)).x;

	bool CirrusDetail = false;

	if (CirrusDetail) {

		CirrusDensity += BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*3.2f,SamplePosition.y)).x;
		CirrusDensity += BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*1.4f,SamplePosition.y)).x;
		CirrusDensity += BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*1.7f,SamplePosition.y*1.2)).x;
		CirrusDensity += BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*1.7f,SamplePosition.y*1.9)).x;
		CirrusDensity += BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*1.7f,SamplePosition.y*2.9)).x;
		CirrusDensity += BetterTexture(u_CirrusClouds, vec2(SamplePosition.x*3.7f,SamplePosition.y*2.9)).x;
		CirrusDensity /= 7.0f;

	}

	const float WeatherWeight = 1.0f;

	CirrusDensity = 1.0f - pow(CirrusDensity, 1.0 / mix(1.0f,15.0f,clamp(u_CirrusStrength * 1.25f,0.0f,1.0f)));
	return (CirrusDensity * 1.0f * 0.2f);
}

float ResolveScatter(float x, float coeff)
{
    float a = -coeff * (1.0 / log(2.0));
    float b = -1.0 / coeff;
    float c =  1.0 / coeff;
    return exp2(a * x) * b + c;
}

// cirrus powder ->
float PlanarPowder(float od)
{
	return 1.0f - exp(-2.0f * od);
}

float CalculatePlanarLightDensity(vec3 P, vec3 D, bool Conservative) {
	
	float RayLength = 15.0f;
    float Accumulated = 0.0;
	int Steps = Conservative ? 4 : 8;

    for(int i = 0; i < Steps; i++)
	{
        Accumulated += RayLength * GetCirrusDensityAt(P);
        P += D * RayLength;
    }

    return Accumulated;

}



// Integrates lighting for cirrus clouds ->

vec4 MarchCirrus(vec3 P, vec3 V, vec2 CosTheta, vec3 SkyColor, vec3 SunLightColor) 
{
	const vec3 Up = normalize(vec3(0.0f,1.0f,0.0f));

	// Base shape ->
	//float BaseOpticalDepth = GetCirrusDensityAt(P) * 1.5f;

	// Base shape ->

	float BaseOpticalDepth = 0.0f;
	vec3 ODPosition = P + V * 0.02f;

	int ShapeSteps = 6;

	for (int Step = 0 ; Step < ShapeSteps ; Step++) {
		
		BaseOpticalDepth += GetCirrusDensityAt(ODPosition) * 1.5f;
		ODPosition += V * 12.0f;
	}

	BaseOpticalDepth /= float(ShapeSteps);

	// Compute transmittance ->
	float Transmittance = exp2(-BaseOpticalDepth); 

	float LightMarch = CalculatePlanarLightDensity(P, NormalizedSUNDIR, false);
	float SkyMarch = CalculatePlanarLightDensity(P, Up, true);

	// Powder ->
	float PowderTerm = PlanarPowder(LightMarch);
    PowderTerm = mix(PowderTerm, 4.0f, MiePhase(CosTheta.x, 0.76f));
    PowderTerm = mix(PowderTerm, 1.0f, CosTheta.x * 0.5f + 0.5f) * PowderTerm;

	float AmbientPowderTerm = PlanarPowder(SkyMarch);
    AmbientPowderTerm = mix(AmbientPowderTerm, 1.0f, CosTheta.y * 0.5f + 0.5f) * AmbientPowderTerm;

	// Phase terms ->
	float DirectPhase = hg(CosTheta.x, LightMarch) * 1.1f;
    float AmbientPhase = hg(CosTheta.y, SkyMarch) * 1.9f;

	// Transmittance terms ->
	float DirectTransmittance = Transmittance * exp(-LightMarch);
	float IndirectTransmittance = Transmittance * exp(-SkyMarch);

	// Combine ->

	float DirectLighting = DirectTransmittance * DirectPhase * PowderTerm * 1.4f;
	float IndirectLighting = IndirectTransmittance * AmbientPhase * AmbientPowderTerm * 1.4f;

	vec3 Scattering = vec3(0.0f);
	Scattering += DirectLighting * (SunLightColor * 0.1f);
	Scattering += IndirectLighting * (SkyColor);

	return vec4(Scattering, Transmittance);
}

const float IsotropicPhaseFunction = 0.25f / PI;
const bool DoFade = true;



void main()
{
	g_TexCoords = v_TexCoords;
	g_TexCoords += (u_JitterValue*3.0f) / u_Dimensions;
	
	
	
	NormalizedSUNDIR = normalize(u_SunDirection);
	o_Data = vec4(0.0f, 0.0f, 0.0f, 1.0f / 2.0f);
	
	vec3 CameraPosition = u_InverseView[3].xyz;
	vec3 Direction = normalize(ComputeRayDirection());

	vec3 AtmosphereAtViewRay = texture(u_Atmosphere, Direction).xyz;

	// We dont need to ray cast the clouds if there is a hit at the current position
	if (!SampleValid(g_TexCoords))
	{
		o_Data.xyz = AtmosphereAtViewRay;
		o_Data.w = 0.2f;
		return;
	}

	int Frame = u_CurrentFrame % 30;
	vec2 BayerIncrement = vec2(Frame * 1.0f, Frame * 0.5f); 
	g_BayerIncrement = BayerIncrement;

	// March from aerian perspective
	vec3 AerialPerspective = vec3(0.0, PlanetRadius + CameraPosition.y, 0.0);

	// Planet sphere ->
	vec2 EarthSphere = RSI(AerialPerspective, Direction, PlanetRadius * 0.99995f);

	// Atmosphere spheres ->
    vec2 SphereMin = RSI(AerialPerspective, Direction, PlanetRadius + CloudAltitudeMin);
    vec2 SphereMax = RSI(AerialPerspective, Direction, PlanetRadius + CloudAltitudeMax);

	if (SphereMin.y >= 49000.0f && DoFade)
	{
		o_Data.xyz = AtmosphereAtViewRay;
		o_Data.w = 0.2f;
		return;
	}

	// Calculate traversal distances ->
	vec2 TraversalDistance = vec2(0.0f);
	TraversalDistance.x = CameraPosition.y > CloudAltitudeMax ? SphereMax.x : SphereMin.y;
	TraversalDistance.y = CameraPosition.y > CloudAltitudeMax ? SphereMin.x : SphereMax.y;

	// March range ->
	float Range = (1.0 - clamp01((CameraPosition.y - CloudAltitudeMax) * 0.1)) * (1.0 - saturate((CloudAltitudeMin - CameraPosition.y) * 0.1));
	
	// Project rays ->
	vec3 ProjectedInner = Direction * TraversalDistance.x;
    ProjectedInner = ProjectedInner * (1.0 - Range);
	vec3 ProjectedOuter = Direction * TraversalDistance.y;
    ProjectedOuter = mix(ProjectedOuter, Direction * CloudAltitudeMin * 12.0, Range);

	// Low frequency hash ->
	float Hash = Bayer32(gl_FragCoord.xy + vec2(u_CurrentFrame * (0.9f/1.0f), u_CurrentFrame * 0.5));
	
	
	int StepCount = clamp(int(u_StepCounts.x),2,64);



	if (CHECKER_STEP_COUNT) {

		bool IsCheckerStep = int(gl_FragCoord.x + gl_FragCoord.y) % 2 == (u_CurrentFrame%2);
		int HalfStepCount = clamp(int(u_StepCounts.x/2),2,64);
		StepCount = clamp(int(mix(HalfStepCount, StepCount, float(IsCheckerStep))),2,64);

	}

	vec3 Increment = (ProjectedOuter - ProjectedInner) / StepCount;
	float StepSize = 1.0f * length(Increment);

	const vec3 AmbientDirection = normalize(vec3(0.0f, 1.0f, 0.0f));
	const vec3 SecondaryBounceDirection = normalize(vec3(0.0f, -1.0f, 0.0f));

	// For phase and powder terms ->
	float VDotL = dot(Direction, NormalizedSUNDIR);
	float VDotAmbient = dot(Direction, AmbientDirection);
	float LDotSecondaryBounce = dot(NormalizedSUNDIR, SecondaryBounceDirection);
	vec3 CosTheta = vec3(VDotL, VDotAmbient, LDotSecondaryBounce);

	// Transmittance, Scattering 
	float Transmittance = 1.0f;

	// x -> DirectScattering, y -> IndirectScattering, z -> Bounced direct scattering 
	vec3 IntegratedScattering = vec3(0.0f);

	// Start from projected inner sphere ->
	vec3 RayPosition = CameraPosition + ProjectedInner + (Increment * Hash);
	RayPosition += Increment * 0.1f;
	
	for (int i = 0 ; i < StepCount ; i++)
	{
	
		if (Transmittance < 0.02f) {
			break;
		}

		float DensityAtPosition = SampleDensity(RayPosition, 0.0f) * 1.2f;
		float CurrentTransmittance = exp2(-DensityAtPosition);
		//float CurrentTransmittance = exp(-DensityAtPosition);

		// Ray march ->
		float DirectDensity = RaymarchLight(RayPosition) * 1.28f;
		float IndirectDensity = RaymarchAmbient(RayPosition) * 1.1f;
		float IndirectSecondBounceDensity = RaymarchBounced(RayPosition) / 2.0f; 

		vec3 ScatterAtPoint = ComputeScattering(RayPosition, DensityAtPosition, CurrentTransmittance, Transmittance, CosTheta, vec3(DirectDensity, IndirectDensity, IndirectSecondBounceDensity));

		IntegratedScattering += ScatterAtPoint;
		Transmittance *= CurrentTransmittance;
		RayPosition += Increment;
	}

	vec3 SunColor = mix(vec3(38.0f),vec3(20.0f),float(1.0f-clamp(u_SunVisibility,0.0f,1.0f)));
	vec3 CirrusSunColor = mix(vec3(38.0f),vec3(40.0f),float(1.0f-clamp(u_SunVisibility,0.0f,1.0f)));
	vec3 SkyColor = pow(texture(u_Atmosphere, normalize(vec3(0.0f, 1.0f, 0.0f))).xyz,vec3(1.0f)) * 4.0f;

	// Add scattering components ->

	vec3 TotalScattering = vec3(0.0f);
	TotalScattering += IntegratedScattering.x * SunColor;
	TotalScattering += IntegratedScattering.y * SkyColor;
	TotalScattering += IntegratedScattering.z * SunColor * clamp(CosTheta.z,0.0f,1.0f) * (1.0f / PI) * 0.5f;
	
	// Cirrus..? 

	float FinalTransmittance = 1.0f;

	FinalTransmittance *= Transmittance;

	// Cirrus -> 

	if (u_CirrusStrength > 0.012f) {
			
		vec3 CirrusPosition = CameraPosition + Direction * SphereMin.y;

		vec4 CirrusMarch = MarchCirrus(CirrusPosition, Direction, vec2(CosTheta.x, CosTheta.y), CirrusSunColor * 0.225f, SkyColor / 3.2525f);

		TotalScattering += CirrusMarch.xyz * FinalTransmittance;

		FinalTransmittance *= CirrusMarch.w;
	}
	

	// Store ->

	vec4 FinalData = vec4(TotalScattering, FinalTransmittance);
	
	if (DoFade) 
	{
		float Fade = 1.0f-(exp(-(SphereMin.y/11000.0f)));
		Fade = clamp(Fade, 0.0f, 1.0f);
		Fade = Fade * Fade;
		Fade = clamp(Fade, 0.0f, 1.0f);
		vec4 Identity = vec4(AtmosphereAtViewRay, 0.3f);
		FinalData = mix(FinalData, Identity, Fade);
	}

	// Output :
	o_Data = vec4(FinalData);

}