#version 330 core

#define SATURATION 1.0f
#define VIBRANCE 1.6f

#define sqr(x) (x*x)
#define square(x) sqr(x)
#define Pow2(x) square(x)
#define Clamp01(x) (clamp(x, 0.0f, 1.0f))
#define clamp01(x) (clamp(x, 0.0f, 1.0f))

//////


// Bayer matrix, used for testing dithering
#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

layout(location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;
in vec3 v_BloomCenter;
flat in int v_PlayerShadowed;

uniform vec3 u_SunDirection;
uniform vec3 u_StrongerLightDirection;
uniform vec2 u_Dimensions;
uniform vec3 u_ViewerPosition;

uniform bool u_SunIsStronger;

uniform bool u_LensFlare = true;
uniform bool u_GodRays = false;
uniform bool u_SSGodRays = false;
uniform bool u_SSAO = false;
uniform bool u_Bloom = false;
uniform bool u_RTAO = false;
uniform bool u_ExponentialFog = false;
uniform bool u_PointVolumetricsToggled = false;


uniform float u_PurkingeEffectStrength;

uniform sampler2D u_Clouds;

uniform bool u_Fucktard;

uniform float u_Time;
uniform float u_FilmGrainStrength;

uniform float u_ChromaticAberrationStrength = 0.0f;

uniform bool u_FilmGrain;

uniform int u_GodRaysStepCount = 12;


uniform int u_CurrentFrame;

uniform bool u_HejlBurgess;

uniform sampler2D u_FramebufferTexture;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;

uniform sampler2D u_VolumetricTexture;
uniform sampler2D u_RTAOTexture;

uniform sampler2D u_BlueNoise;
uniform sampler2D u_SSAOTexture;

uniform sampler2D u_BloomMips[5];
uniform sampler2D u_BloomBrightTexture;

uniform sampler2D u_ShadowTexture;


uniform sampler2D u_CloudData;


uniform sampler2D u_VolumetricsCompute;
uniform sampler2D u_LensDirtOverlay;
uniform samplerCube u_Sky;

uniform samplerCube u_NightSkymap;


uniform mat4 u_ProjectionMatrix;
uniform mat4 u_ViewMatrix;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;


uniform mat4 u_RotationMatrix;

uniform float u_LensFlareIntensity;
uniform float u_Exposure;

uniform float u_NebulaStrength;

uniform float u_GodRaysStrength;

uniform float u_BloomStrength;

uniform bool u_LensDirt;
uniform float u_LensDirtStrength;
uniform bool u_HQBloomUpscale;

vec4 textureBicubic(sampler2D sampler, vec2 texCoords);


vec4 TextureSmooth(sampler2D samp, vec2 uv) 
{
    vec2 textureResolution = textureSize(samp, 0).xy;
	uv = uv*textureResolution + 0.5f;
	vec2 iuv = floor(uv);
	vec2 fuv = fract(uv);
	uv = iuv + fuv*fuv*(3.0f-2.0f*fuv); 
	uv = (uv - 0.5f) / textureResolution;
	return texture(samp, uv).xyzw;
}

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SamplePositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

bool CompareFloatNormal(float x, float y) {
    return abs(x - y) < 0.02f;
}

vec3 GetNormalFromID(float n) {
	const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f));
    int idx = int(round(n*10.0f));

    if (idx > 5) {
        return vec3(1.0f, 1.0f, 1.0f);
    }

    return Normals[idx];
}

vec3 SampleNormal(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}

void Vignette(inout vec3 color) 
{
	float dist = distance(v_TexCoords.xy, vec2(0.5f)) * 2.0f;
	dist /= 1.5142f;
	color.rgb *= 1.0f - dist * 0.5;
}

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

vec2 WorldToScreen(vec3 pos)
{
    vec4 ViewSpace = u_ViewMatrix * vec4(pos, 1.0f);
    vec4 Projected = u_ProjectionMatrix * ViewSpace;
    Projected.xyz /= Projected.w;
    Projected.xyz = Projected.xyz * 0.5f + 0.5f;

    return Projected.xy;
} 

float GetScreenSpaceGodRays(vec3 position)
{
	vec3 view_position = (u_ViewMatrix * vec4(position, 1.0f)).xyz;

    vec2 SunScreenSpacePosition = WorldToScreen(u_StrongerLightDirection * 10000.0f); 

    float ScreenSpaceDistToSun = length(v_TexCoords - SunScreenSpacePosition.xy);
    float RayIntensity = clamp(1.0f - ScreenSpaceDistToSun, 0.0f, 0.75f);
    float RayIntensityMultiplier = u_StrongerLightDirection == u_SunDirection ? 0.35224f : 0.15f;

    float rays = 0.0;
    int SAMPLES = clamp(u_GodRaysStepCount, 8, 65);
	float dither = texture(u_BlueNoise, v_TexCoords * (u_Dimensions / vec2(256.0f))).r;

    for (int i = 0; i < SAMPLES; i++)
    {
        float scale = (1.0f - (float(i) / float(SAMPLES))) + dither / float(SAMPLES);

        vec2 coord = (v_TexCoords - SunScreenSpacePosition) * scale + SunScreenSpacePosition;
        coord = clamp(coord, 0.001f, 0.999f);

        float is_sky_at = SamplePositionAt(u_PositionTexture, coord).w > 0.0f ? 0.0f : 1.0f;

        rays += is_sky_at * RayIntensity * RayIntensityMultiplier;
    }

	rays /=  float(SAMPLES);

    return rays;
}

// not my code.
// by mu6k

float Fnoise(float t)
{
	return texture(u_BlueNoise, vec2(t, 0.0f) / vec2(256).xy).x;
}

float Fnoise(vec2 t)
{
	return texture(u_BlueNoise, t / vec2(256).xy).x;
}

vec3 lensflare(vec2 uv, vec2 pos)
{
	vec2 main = uv - pos;
	vec2 uvd = uv * (length(uv));
	float ang = atan(main.x,main.y);
	float dist = length(main); 
	dist = pow(dist, 0.1f);
	float n = Fnoise(vec2(ang * 16.0f, dist * 32.0f));
	float f0 = 1.0 / (length(uv - pos) * 16.0f + 1.0f);
	f0 = f0 + f0 * (sin(Fnoise(sin(ang * 2.0f + pos.x) * 4.0f - cos(ang * 3.0f + pos.y)) * 16.0f) * 0.1 + dist * 0.1f + 0.8f);
	float f1 = max(0.01f - pow(length(uv + 1.2f * pos), 1.9f), 0.0f) * 7.0f;
	float f2 = max(1.0f / (1.0f + 32.0f * pow(length(uvd + 0.8f * pos), 2.0f)), 0.0f) * 0.25f;
	float f22 = max(1.0f / (1.0f + 32.0f * pow(length(uvd + 0.85f * pos), 2.0f)), 0.0f) * 0.23f;
	float f23 = max(1.0f / (1.0f + 32.0f * pow(length(uvd + 0.9f * pos), 2.0f)), 0.0f) * 0.21f;
	vec2 uvx = mix(uv, uvd, -0.5f);
	float f4 = max(0.01f - pow(length(uvx + 0.4f * pos), 2.4f), 0.0f) * 6.0f;
	float f42 = max(0.01f - pow(length(uvx + 0.45f * pos), 2.4f), 0.0f) * 5.0f;
	float f43 = max(0.01f - pow(length(uvx + 0.5f *pos), 2.4f), 0.0f) * 3.0f;
	uvx = mix(uv, uvd, -0.4f);
	float f5 = max(0.01f - pow(length(uvx + 0.2f * pos), 5.5f), 0.0f) * 2.0f;
	float f52 = max(0.01f - pow(length(uvx + 0.4f * pos), 5.5f), 0.0f) * 2.0f;
	float f53 = max(0.01f - pow(length(uvx + 0.6f * pos), 5.5f), 0.0f) * 2.0f;
	uvx = mix(uv, uvd, -0.5f);
	float f6 = max(0.01f - pow(length(uvx - 0.3f * pos), 1.6f), 0.0f) * 6.0f;
	float f62 = max(0.01f - pow(length(uvx - 0.325f * pos), 1.6f), 0.0f) * 3.0f;
	float f63 = max(0.01f - pow(length(uvx - 0.35f * pos), 1.6f), 0.0f) * 5.0f;
	vec3 c = vec3(0.0f);
	c.r += f2 + f4 + f5 + f6; 
	c.g += f22 + f42 + f52 + f62;
	c.b += f23 + f43 + f53 + f63;
	c = c * 1.3f - vec3(length(uvd) * 0.05f);
	return c;
}

vec3 ToClipSpace(in vec3 pos)
{
	vec4 ViewSpace = u_ViewMatrix * vec4(pos, 1.0f);
    vec4 Projected = u_ProjectionMatrix * ViewSpace;
    Projected.xyz /= Projected.w;

    return Projected.xyz;
}

// Exponential Distance fog
vec3 ApplyFog(in vec3 InputColor, in float dist_to_point)
{
	dist_to_point *= 2.0f;
	float b = 0.00875f;
    float FogAmount = 1.0 - exp(-dist_to_point * b);
    vec3 FogColor  = texture(u_Sky, normalize(v_RayDirection)).xyz * 1.5f; //vec3(0.5f, 0.6f, 0.7f);
    return mix(InputColor, FogColor, FogAmount );
}

// Exponential height fog
vec3 ApplyFog(in vec3 InputColor, in float dist_to_point, in vec3 RayDir, in vec3 SunDir)  
{
	float b = 0.00575f;
    float FogAmount = 1.0 - exp(-dist_to_point * b);
    float SunAmount = max(dot(RayDir, SunDir), 0.0);
    vec3  FogColor  = mix(vec3(0.5f, 0.6f, 0.7f), vec3(1.0f, 0.9f, 0.7f), pow(SunAmount , 8.0f));
    return mix(InputColor, FogColor, FogAmount );
}


float Brightness(vec3 c)
{
    return max(max(c.r, c.g), c.b);
}



// Chromatic aberration.
vec3 BasicChromaticAberation()
{
	const float ChannelCount = 3.0f;
	float AberationScale = mix(0.0f, 0.5f, u_ChromaticAberrationStrength);
	vec2 DistanceWeight = v_TexCoords - 0.5f;
    vec2 Aberrated = AberationScale * pow(DistanceWeight, vec2(ChannelCount));
    vec3 Final = vec3(0.0f);
	float TotalWeight = 0.01f;
    
    for (int i = 1; i <= 12; i++)
    {
        float wg = 1.0f / pow(2.0f, float(i)); // Blur Weights, tested.

		if (v_TexCoords - float(i) * Aberrated == clamp(v_TexCoords - float(i) * Aberrated, 0.0001f, 0.9999f)) {
			Final.r += texture(u_FramebufferTexture, v_TexCoords - float(i) * Aberrated).r * wg;
		}

		if (v_TexCoords + float(i) * Aberrated == clamp(v_TexCoords + float(i) * Aberrated, 0.0001f, 0.9999f)) {
			Final.b += texture(u_FramebufferTexture, v_TexCoords + float(i) * Aberrated).b * wg;
		}

		TotalWeight += wg;
    }
    
	//const float TotalWeight = 0.9961f; //(1.0 / pow(2.0f, float(i)) i = 1 -> 8 
	Final.g = texture(u_FramebufferTexture, v_TexCoords).g * TotalWeight;
	return max(Final,0.0f);
}

const mat3x3 xyzToRGBMatrix = mat3(
    3.1338561, -1.6168667, -0.4906146,
    -0.9787684,  1.9161415,  0.0334540,
    0.0719453, -0.2289914,  1.4052427
);

const mat3x3 rgbToXYZMatrix = mat3(
    vec3(0.5149, 0.3244, 0.1607),
    vec3(0.3654, 0.6704, 0.0642),
    vec3(0.0248, 0.1248, 0.8504)
);

vec3 xyzToRGB(in vec3 xyz) {
    float r = dot(xyz, xyzToRGBMatrix[0]);
    float g = dot(xyz, xyzToRGBMatrix[1]);
    float b = dot(xyz, xyzToRGBMatrix[2]);
    return vec3(r, g, b);
}

vec3 rgbToXYZ(in vec3 rgb) {
    float x = dot(rgb, rgbToXYZMatrix[0]);
    float y = dot(rgb, rgbToXYZMatrix[1]);
    float z = dot(rgb, rgbToXYZMatrix[2]);
    return vec3(x, y, z);
}

vec3 PurkinjeEffect(vec3 Color) 
{
	vec3 BaseColor = Color;
    vec3 RodResponse = vec3(7.15e-5f, 4.81e-1f, 3.28e-1f);
    vec3 XYZ = rgbToXYZ(Color);
    vec3 ScopticLuma = XYZ * (1.33f * (1.0f + (XYZ.y + XYZ.z) / XYZ.x) - 1.68f);
    float Purkinge = dot(RodResponse, xyzToRGB(ScopticLuma));
    Color = mix(Color, Purkinge * vec3(0.5f, 0.7f, 1.0f), exp2(-Purkinge * 2.0f));
    return mix(max(Color, 0.0), BaseColor, clamp(1.0f-u_PurkingeEffectStrength,0.0f,1.0f));
}



// Used to get a color for a temperature ->
vec3 Blackbody(float temperature) 
{ 
	const mat3 tmp_r2020 = mat3(
		0.708 / 0.292, 0.17  / 0.797, 0.131 / 0.046,
		1.0,           1.0,           1.0,
		0.0   / 0.292, 0.033 / 0.797, 0.823 / 0.046
	);
	const vec3 lc_r2020 = vec3(0.3127 / 0.3290, 1.0, 0.3583 / 0.3290) * inverse(tmp_r2020);
	const mat3 R2020ToXyz = mat3(lc_r2020 * tmp_r2020[0], lc_r2020, lc_r2020 * tmp_r2020[2]);
	const mat3 XyzToR2020 = inverse(R2020ToXyz);
	const mat3 XyzToRgb = XyzToR2020;

	const vec4[2] xc = vec4[2](
		vec4(-0.2661293e9,-0.2343589e6, 0.8776956e3, 0.179910), // 1667k <= t <= 4000k
		vec4(-3.0258469e9, 2.1070479e6, 0.2226347e3, 0.240390)  // 4000k <= t <= 25000k
	);

	const vec4[3] yc = vec4[3](
		vec4(-1.1063814,-1.34811020, 2.18555832,-0.20219683), // 1667k <= t <= 2222k
		vec4(-0.9549476,-1.37418593, 2.09137015,-0.16748867), // 2222k <= t <= 4000k
		vec4( 3.0817580,-5.87338670, 3.75112997,-0.37001483)  // 4000k <= t <= 25000k
	);

	float temperatureSquared = temperature * temperature;
	vec4 t = vec4(temperatureSquared * temperature, temperatureSquared, temperature, 1.0);
	float x = dot(1.0 / t, temperature < 4000.0 ? xc[0] : xc[1]);
	float xSquared = x * x;
	vec4 xVals = vec4(xSquared * x, xSquared, x, 1.0);
	vec3 xyz = vec3(0.0);
	xyz.y = 1.0;
	xyz.z = 1.0 / dot(xVals, temperature < 2222.0 ? yc[0] : temperature < 4000.0 ? yc[1] : yc[2]);
	xyz.x = x * xyz.z;
	xyz.z = xyz.z - xyz.x - 1.0;

	return xyz * XyzToRgb;
}

float Noise2d(in vec2 x)
{
    float xhash = cos(x.x * 37.0);
    float yhash = cos(x.y * 57.0);
    return fract(415.92653 * (xhash + yhash));
}

// thresholded white noise :
float NoisyStarField(in vec2 UV, float Threshold)
{
    float StarVal = Noise2d(UV);

    if (StarVal >= Threshold) 
    {
        StarVal = pow((StarVal - Threshold) / (1.0f - Threshold), 6.0f);
    }

    else 
    {
        StarVal = 0.0;
    }

    return StarVal;
}

vec2 Hash2(vec3 p3) {
	p3 = fract(p3 * vec3(443.897, 441.423, 437.195));
	p3 += dot(p3, p3.yzx + 19.19);
	return fract((p3.xx + p3.yz) * p3.zy);
}

float LinearStep(float e0, float e1, float x) { return clamp01((x - e0) / (e1 - e0)); }

vec2 ProjectDirection(vec3 Direction, vec2);

vec3 ShadeStars(vec3 sky, vec3 Lo, float sv, float base_transmittance)
{
    const float StarScale = 256.0f;
	const float FakeCoverage = 0.01f;
	const float LuminanceMax = 0.7f * 5.0f;
	const float BlackbodyMinTemperature = 1000.0;
	const float BlackbodyMaxTemperature = 9000.0;

	vec3  p = Lo * StarScale;
	ivec3 i = ivec3(floor(p));
	vec3  f = p - i;
	float r = dot(f - 0.5f, f - 0.5f);

	vec2 hash = Hash2(i);
	hash.y = 2.0f * hash.y - 4.0f * hash.y * hash.y + 3.0f * hash.y * hash.y * hash.y;

	vec3 BlackbodyColor = pow(Blackbody(mix(BlackbodyMinTemperature, BlackbodyMaxTemperature, hash.y)), vec3(1.33f));
	vec3 Stars = clamp((LuminanceMax * LinearStep(0.25f, 0.0f, r) * sqr(LinearStep(1.0f-FakeCoverage, 1.0f, hash.x)) * BlackbodyColor * sv * base_transmittance), 0.0f, 8.0f);
	float StarTransmittance = 1.0f; // GetLuminance(Stars) * LinearStep(0.2,0.6f,GetLuminance(BlackbodyColor)) * sqr(base_transmittance) * sv 
	return sky * StarTransmittance + Stars;
}



vec4 SampleTextureCatmullRom(sampler2D tex, in vec2 uv); // catmull rom texture interp

 #define VOLUMETRIC_BICUBIC

vec4 texture_catmullrom(sampler2D tex, vec2 uv);
void main()
{

	vec4 PositionAt = SamplePositionAt(u_PositionTexture, v_TexCoords).rgba;
	vec3 NormalAt = SampleNormal(u_NormalTexture, v_TexCoords).rgb;
	vec3 ClipSpaceAt = ToClipSpace(PositionAt.xyz);
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	float BayerDither = Bayer64(gl_FragCoord.xy);
	vec2 InverseVolRes = 1.0f / textureSize(u_VolumetricsCompute,0);
	vec3 PointVolumetrics = vec3(0.0f);
	
	if (u_PointVolumetricsToggled) {
	
		PointVolumetrics.xyz = textureBicubic(u_VolumetricsCompute, v_TexCoords+(Bayer128(gl_FragCoord.xy)*1.2*(1.0f/textureSize(u_VolumetricsCompute,0)))).xyz;
	}


	if ( (PositionAt.w > 0.0f ))
	{
		vec3 InputColor;
		InputColor = u_ChromaticAberrationStrength <= 0.001f ? texture(u_FramebufferTexture, v_TexCoords).rgb : BasicChromaticAberation();

		if (u_SSAO && (!u_RTAO))
		{
			// Use non linear projected depth to fade ssao :
			float ssao_strength = 0.0f;
			float max_ssao_strength = 2.5f;
			ssao_strength = ((1.0f - ClipSpaceAt.z) * 100.0f) * max_ssao_strength;
			ssao_strength = clamp(ssao_strength, 0.0f, max_ssao_strength);
			//float SampledSSAO = DepthOnlyBilateralUpsample(u_SSAOTexture, v_TexCoords, PositionAt.z).r;
			float SampledSSAO = texture(u_SSAOTexture, v_TexCoords).r;
			float SSAO = pow(SampledSSAO, ssao_strength);
			SSAO = clamp(SSAO, 0.00001, 1.0f);
			InputColor *= SSAO;
		}

		if (u_RTAO)
		{
			float rtao_strength = 0.0f;
			float max_rtao_strength = 1.1f;
			rtao_strength = ((1.0f - ClipSpaceAt.z) * 200.0f) * max_rtao_strength;
			rtao_strength = clamp(rtao_strength, 0.0f, max_rtao_strength);
			//float RTAO = BilateralUpsample(u_RTAOTexture, v_TexCoords, NormalAt, PositionAt.w).r;
			float RTAO = texture(u_RTAOTexture, v_TexCoords).r;
			RTAO = pow(RTAO, rtao_strength);
			RTAO = max(RTAO, 0.825f);
			InputColor = (RTAO * InputColor) + (0.01f * InputColor);
		}



		if (u_GodRays)
		{
			float fake_vol_multiplier = u_SunIsStronger ? 1.21f : 0.9f;
			float intensity = u_SunIsStronger ? 0.55f : 0.025f;
			//float god_rays = DepthOnlyBilateralUpsample(u_VolumetricTexture, v_TexCoords, PositionAt.z).r * intensity * u_GodRaysStrength;
			float god_rays = texture(u_VolumetricTexture, v_TexCoords).r * intensity * u_GodRaysStrength;
			vec3 ss_volumetric_color = u_SunIsStronger ? (vec3(250.0f, 250.0f, 230.0f) / 255.0f) : (vec3(96.0f, 192.0f, 255.0f) / 255.0f);
			InputColor += god_rays * ss_volumetric_color;
			fake_vol_multiplier = 0.3f;
		}

		if (u_SSGodRays)
		{
			float fake_vol_multiplier = u_SunIsStronger ? 1.21f : 0.9f;
			float god_rays = GetScreenSpaceGodRays(PositionAt.xyz) * fake_vol_multiplier * u_GodRaysStrength;
			vec3 ss_volumetric_color = u_SunIsStronger ? (vec3(250.0f, 250.0f, 230.0f) / 255.0f) : (vec3(96.0f, 192.0f, 255.0f) / 255.0f);
			InputColor += god_rays * ss_volumetric_color;
		}

		InputColor += PointVolumetrics;


		if (false)
		{
			InputColor = ApplyFog(InputColor, PositionAt.w); 
			//InputColor = ApplyFog(InputColor, PositionAt.w, v_RayDirection, normalize(u_StrongerLightDirection));
		}

		o_Color = InputColor;

	}

	else 
	{
		vec3 Lo = normalize(v_RayDirection);
		float star_visibility;
		star_visibility = 1.- clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; 
		
		vec3 Nighteffects = vec3(0.0f);
		if (star_visibility > 0.01f&&u_NebulaStrength>0.01f) {
			Nighteffects = texture(u_NightSkymap, (u_RotationMatrix*vec4(Lo.xyz,1.)).xyz*2.).xyz;
		}

		float transmittance = texture(u_Clouds, v_TexCoords).w;


		o_Color = u_ChromaticAberrationStrength <= 0.001f ? texture(u_FramebufferTexture, v_TexCoords).rgb : BasicChromaticAberation() ;
		o_Color = o_Color  + ((Nighteffects * 0.5f * u_NebulaStrength * sqr(star_visibility) * transmittance));
		o_Color = ShadeStars(o_Color, Lo, star_visibility, transmittance * transmittance);
		o_Color += PointVolumetrics;
	}



	
	if (u_Bloom)
	{
		vec3 Bloom[5] = vec3[](vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f));

		// bicubic upsampling because the bloom textures are super low res
		// res : 0.25, 0.125, 0.0625, 0.03125
		// 4 * 4 + 9 (or 1 if blur is disabled) samples 
		// We sample the bright texture to preserve some detail from the emissive textures 

		vec3 BaseBrightTex = vec3(0.0f);

		bool Blur = false;

		if (Blur) {

			vec2 TexelSizeBlur = 1./textureSize(u_BloomBrightTexture,0).xy;
			
			for (int x = -1 ; x <= 1 ; x++) {
				for (int y = -1 ; y <= 1 ; y++) {
					BaseBrightTex += texture(u_BloomBrightTexture, v_TexCoords + vec2(x, y) * TexelSizeBlur * 1.0f).xyz;
				}
			}

			BaseBrightTex *= 1.0f / 9.0f;
		}

		else {
			BaseBrightTex = textureBicubic(u_BloomBrightTexture, v_TexCoords).xyz;
		}



		if (u_HQBloomUpscale) {
			Bloom[0] += textureBicubic(u_BloomMips[0], v_TexCoords).xyz;
			Bloom[1] += textureBicubic(u_BloomMips[1], v_TexCoords).xyz; 
			Bloom[2] += textureBicubic(u_BloomMips[2], v_TexCoords).xyz; 
			Bloom[3] += textureBicubic(u_BloomMips[3], v_TexCoords).xyz; 
			Bloom[4] += textureBicubic(u_BloomMips[4], v_TexCoords).xyz; 
		}

		else {
			Bloom[0] += TextureSmooth(u_BloomMips[0], v_TexCoords).xyz;
			Bloom[1] += TextureSmooth(u_BloomMips[1], v_TexCoords).xyz; 
			Bloom[2] += TextureSmooth(u_BloomMips[2], v_TexCoords).xyz; 
			Bloom[3] += TextureSmooth(u_BloomMips[3], v_TexCoords).xyz; 
			Bloom[4] += TextureSmooth(u_BloomMips[4], v_TexCoords).xyz; 
		}

		vec3 TotalBloom = vec3(0.0f);

		// Tweaks the bloom color a bit
		float AmplificationFactor = 1.1;
			
		// Weighted average ->

		float Weights[5] = float[5](5.75f, 4.95f, 4.9f, 4.8f, 4.75f);
		const float DetailWeight = 7.2f;

		TotalBloom = (BaseBrightTex * DetailWeight * 1.0f) + TotalBloom;
		TotalBloom = (pow(Bloom[0], vec3(1.0f / AmplificationFactor)) * Weights[0]) + TotalBloom;
		TotalBloom = (pow(Bloom[1], vec3(1.0f / AmplificationFactor)) * Weights[1]) + TotalBloom;
		TotalBloom = (pow(Bloom[2], vec3(1.0f / AmplificationFactor)) * Weights[2]) + TotalBloom;
		TotalBloom = (pow(Bloom[3], vec3(1.0f / AmplificationFactor)) * Weights[3]) + TotalBloom;
		TotalBloom = (pow(Bloom[4], vec3(1.0f / AmplificationFactor)) * Weights[4]) + TotalBloom;

		float TotalWeights = DetailWeight + Weights[0] + Weights[1] + Weights[2] + Weights[3] + Weights[4];
		TotalBloom /= TotalWeights;

		if (u_LensDirt) {
			vec3 LensDirtFetch = TextureSmooth(u_LensDirtOverlay, v_TexCoords).xyz;
			TotalBloom += v_BloomCenter * 5.33f * LensDirtFetch * u_LensDirtStrength * pow(1.0f-distance(v_TexCoords, vec2(0.5f)),6.0f);
		}

		o_Color += TotalBloom * u_BloomStrength ;
	}

	else 
	{
		//float Emissivity = texture(u_PBRTexture, v_TexCoords).w;
		//o_Color += (Emissivity * 8.0f) * o_Color;
	}

	if (u_LensFlare && u_SunIsStronger && v_PlayerShadowed==0)
	{
		vec2 SunScreenSpacePosition = WorldToScreen(u_SunDirection * 10000.0f) - 0.5f; 
		SunScreenSpacePosition.x *= u_Dimensions.x / u_Dimensions.y;
		vec2 LensFlareCoord = v_TexCoords - 0.5f;
		LensFlareCoord.x *= u_Dimensions.x / u_Dimensions.y;
		
		vec3 LensFlare = vec3(1.6f, 1.2f, 1.0f) * lensflare(LensFlareCoord, SunScreenSpacePosition);
		LensFlare = clamp(LensFlare, 0.0f, 0.9999f);
		o_Color += LensFlare * u_LensFlareIntensity;
	}




	if (u_PurkingeEffectStrength > 0.01f) {
		o_Color.xyz = PurkinjeEffect(o_Color.xyz);
	}
}


vec4 cubic(float v){
    vec4 n = vec4(1.0, 2.0, 3.0, 4.0) - v;
    vec4 s = n * n * n;
    float x = s.x;
    float y = s.y - 4.0 * s.x;
    float z = s.z - 4.0 * s.y + 6.0 * s.x;
    float w = 6.0 - x - y - z;
    return vec4(x, y, z, w) * (1.0/6.0);
}

vec4 textureBicubic(sampler2D sampler, vec2 texCoords)
{

   vec2 texSize = textureSize(sampler, 0);
   vec2 invTexSize = 1.0 / texSize;

   texCoords = texCoords * texSize - 0.5;


    vec2 fxy = fract(texCoords);
    texCoords -= fxy;

    vec4 xcubic = cubic(fxy.x);
    vec4 ycubic = cubic(fxy.y);

    vec4 c = texCoords.xxyy + vec2 (-0.5, +1.5).xyxy;

    vec4 s = vec4(xcubic.xz + xcubic.yw, ycubic.xz + ycubic.yw);
    vec4 offset = c + vec4 (xcubic.yw, ycubic.yw) / s;

    offset *= invTexSize.xxyy;

    vec4 sample0 = texture(sampler, offset.xz);
    vec4 sample1 = texture(sampler, offset.yz);
    vec4 sample2 = texture(sampler, offset.xw);
    vec4 sample3 = texture(sampler, offset.yw);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix(
       mix(sample3, sample2, sx), mix(sample1, sample0, sx)
    , sy);
}

#define sqr(x) (x*x)

vec2 ClampUV(vec2 x) {
	return clamp(x, 0.000001f, 0.999999f);
}


vec4 texture_catmullrom(sampler2D tex, vec2 uv) {
    vec2 res    = textureSize(tex, 0);
	vec2 pixelSize = 1.0f / res;

    vec2 coord  = uv*res;
    vec2 coord1 = floor(coord - 0.5) + 0.5;

    vec2 f      = coord-coord1;

    vec2 w0     = f*(-0.5 + f*(1.0-0.5*f));
    vec2 w1     = 1.0 + sqr(f)*(-2.5+1.5*f);
    vec2 w2     = f*(0.5 + f*(2.0-1.5*f));
    vec2 w3     = sqr(f)*(-0.5+0.5*f);

    vec2 w12    = w1+w2;
    vec2 delta12 = w2/w12;

    vec2 uv0    = coord1 - vec2(1.0);
    vec2 uv3    = coord1 + vec2(1.0);
    vec2 uv12   = coord1 + delta12;

        uv0    *= pixelSize;
        uv3    *= pixelSize;
        uv12   *= pixelSize;

    vec4 col    = vec4(0.0);
        col    += textureLod(tex, ClampUV(vec2(uv0.x, uv0.y)), 0)*w0.x*w0.y;
        col    += textureLod(tex, ClampUV(vec2(uv12.x, uv0.y)), 0)*w12.x*w0.y;
        col    += textureLod(tex, ClampUV(vec2(uv3.x, uv0.y)), 0)*w3.x*w0.y;
        col    += textureLod(tex, ClampUV(vec2(uv0.x, uv12.y)), 0)*w0.x*w12.y;
        col    += textureLod(tex, ClampUV(vec2(uv12.x, uv12.y)), 0)*w12.x*w12.y;
        col    += textureLod(tex, ClampUV(vec2(uv3.x, uv12.y)), 0)*w3.x*w12.y;
        col    += textureLod(tex, ClampUV(vec2(uv0.x, uv3.y)), 0)*w0.x*w3.y;
        col    += textureLod(tex, ClampUV(vec2(uv12.x, uv3.y)), 0)*w12.x*w3.y;
        col    += textureLod(tex, ClampUV(vec2(uv3.x, uv3.y)), 0)*w3.x*w3.y;

    return clamp(col, 0.0, 65535.0);
}


vec2 ProjectDirection(vec3 Direction, vec2 TextureSize) 
{
    //vec2 TextureSize = vec2(textureSize(u_DebugTexture,0).xy);
	float TileSize = min(floor(TextureSize.x * 0.5f) / 1.5f, floor(TextureSize.y * 0.5f));
	float TileSizeDivided = (0.5f * TileSize) - 1.5f;
	vec2 CurrentCoordinate;

	if (abs(Direction.x) > abs(Direction.y) && abs(Direction.x) > abs(Direction.z)) 
    {
		Direction /= max(abs(Direction.x), 0.001f);
		CurrentCoordinate.x = Direction.y * TileSizeDivided + TileSize * 0.5f;
		CurrentCoordinate.y = Direction.z * TileSizeDivided + TileSize * (Direction.x < 0.0f ? 0.5f : 1.5f);
	} 
    
    else if (abs(Direction.y) > abs(Direction.x) && abs(Direction.y) > abs(Direction.z))
    {
		Direction /= max(abs(Direction.y), 0.001f);
		CurrentCoordinate.x = Direction.x * TileSizeDivided + TileSize * 1.5f;
		CurrentCoordinate.y = Direction.z * TileSizeDivided + TileSize * (Direction.y < 0.0f ? 0.5f : 1.5f);
	} 
    
    else 
    {
		Direction /= max(abs(Direction.z), 0.001f);
		CurrentCoordinate.x = Direction.x * TileSizeDivided + TileSize * 2.5f;
		CurrentCoordinate.y = Direction.y * TileSizeDivided + TileSize * (Direction.z < 0.0f ? 0.5f : 1.5f);
	}

	return CurrentCoordinate / max(TextureSize, 0.01f);
}