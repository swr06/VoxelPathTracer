#version 330 core




//////
#define CG_RR 255 
#define CG_RG 0 
#define CG_RB 0 
#define CG_RI 1.00 
#define CG_RM 0 
#define CG_RC 1.00 

#define CG_GR 0 
#define CG_GG 255 
#define CG_GB 0 
#define CG_GI 1.00 
#define CG_GM 0 
#define CG_GC 1.00 

#define CG_BR 0 
#define CG_BG 0 
#define CG_BB 255
#define CG_BI 1.00 
#define CG_BM 0 
#define CG_BC 1.00 

#define CG_TR 255 
#define CG_TG 255 
#define CG_TB 255 
#define CG_TI 1.00 
#define CG_TM 0.0 

#define SATURATION 1.0f
#define VIBRANCE 1.6f
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

uniform sampler2D u_BloomMips[4];

uniform sampler2D u_ShadowTexture;

uniform sampler2D u_PBRTexture;

uniform sampler2D u_CloudData;


uniform sampler2D u_VolumetricsCompute;

uniform mat4 u_ProjectionMatrix;
uniform mat4 u_ViewMatrix;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform float u_LensFlareIntensity;
uniform float u_Exposure;

uniform float u_GodRaysStrength;

vec4 textureBicubic(sampler2D sampler, vec2 texCoords);

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

//Due to low sample count we "tonemap" the inputs to preserve colors and smoother edges
vec3 WeightedSample(sampler2D colorTex, vec2 texcoord)
{
	vec3 wsample = texture(colorTex,texcoord).rgb * 1.0f;
	return wsample / (1.0f + GetLuminance(wsample));
}

vec3 smoothfilter(in sampler2D tex, in vec2 uv)
{
	vec2 textureResolution = textureSize(tex, 0);
	uv = uv*textureResolution + 0.5;
	vec2 iuv = floor( uv );
	vec2 fuv = fract( uv );
	uv = iuv + fuv*fuv*fuv*(fuv*(fuv*6.0-15.0)+10.0);
	uv = (uv - 0.5)/textureResolution;
	return WeightedSample( tex, uv);
}

vec3 sharpen(in sampler2D tex, in vec2 coords) 
{
	vec2 renderSize = textureSize(tex, 0);
	float dx = 1.0 / renderSize.x;
	float dy = 1.0 / renderSize.y;
	vec3 sum = vec3(0.0);
	sum += -1. * texture(tex, coords + vec2( -1.0 * dx , 0.0 * dy)).rgb;
	sum += -1. * texture(tex, coords + vec2( 0.0 * dx , -1.0 * dy)).rgb;
	sum += 5. * texture(tex, coords + vec2( 0.0 * dx , 0.0 * dy)).rgb;
	sum += -1. * texture(tex, coords + vec2( 0.0 * dx , 1.0 * dy)).rgb;
	sum += -1. * texture(tex, coords + vec2( 1.0 * dx , 0.0 * dy)).rgb;
	return sum;
}

mat3 ACESInputMat = mat3(
    0.59719, 0.07600, 0.02840,
    0.35458, 0.90834, 0.13383,
    0.04823, 0.01566, 0.83777
);

// ODT_SAT => XYZ => D60_2_D65 => sRGB
mat3 ACESOutputMat = mat3(
    1.60475, -0.10208, -0.00327,
    -0.53108, 1.10813, -0.07276,
    -0.07367, -0.00605, 1.07602
);

vec3 RRTAndODTFit(vec3 v)
{
    vec3 a = v * (v + 0.0245786f) - 0.000090537f;
    vec3 b = v * (0.983729f * v + 0.4329510f) + 0.238081f;
    return a / b;
}

vec4 ACESFitted(vec4 Color, float Exposure)
{
    Color.rgb *= Exposure * 0.6;
    
    Color.rgb = ACESInputMat * Color.rgb;
    Color.rgb = RRTAndODTFit(Color.rgb);
    Color.rgb = ACESOutputMat * Color.rgb;

    return Color;
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

void ColorGrading(inout vec3 color)
{
	vec3 cgColor = pow(color.r, CG_RC) * pow(vec3(CG_RR, CG_RG, CG_RB) / 255.0, vec3(2.2)) +
				   pow(color.g, CG_GC) * pow(vec3(CG_GR, CG_GG, CG_GB) / 255.0, vec3(2.2)) +
				   pow(color.b, CG_BC) * pow(vec3(CG_BR, CG_BG, CG_BB) / 255.0, vec3(2.2));
	vec3 cgMin = pow(vec3(CG_RM, CG_GM, CG_BM) / 255.0, vec3(2.2));
	color = (cgColor * (1.0 - cgMin) + cgMin) * vec3(CG_RI, CG_GI, CG_BI);
	
	vec3 cgTint = pow(vec3(CG_TR, CG_TG, CG_TB) / 255.0, vec3(2.2)) * GetLuminance(color) * CG_TI;
	color = mix(color, cgTint, CG_TM);
}

void ColorSaturation(inout vec3 color) 
{
	float grayVibrance = (color.r + color.g + color.b) / 3.0;
	float graySaturation = grayVibrance;
	if (SATURATION < 1.00) graySaturation = dot(color, vec3(0.299, 0.587, 0.114));

	float mn = min(color.r, min(color.g, color.b));
	float mx = max(color.r, max(color.g, color.b));
	float sat = (1.0 - (mx - mn)) * (1.0 - mx) * grayVibrance * 5.0;
	vec3 lightness = vec3((mn + mx) * 0.5);

	color = mix(color, mix(color, lightness, 1.0 - VIBRANCE), sat);
	color = mix(color, lightness, (1.0 - lightness) * (2.0 - VIBRANCE) / 2.0 * abs(VIBRANCE - 1.0));
	color = color * SATURATION - graySaturation * (SATURATION - 1.0);
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

float Fnoise(float t)
{
	return texture(u_BlueNoise, vec2(t, 0.0f) / vec2(256).xy).x;
}

float Fnoise(vec2 t)
{
	return texture(u_BlueNoise, t / vec2(256).xy).x;
}

// by mu6k
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

vec3 DepthOnlyBilateralUpsample(sampler2D tex, vec2 txc, float base_depth)
{
    const vec2 Kernel[4] = vec2[](
        vec2(0.0f, 1.0f),
        vec2(1.0f, 0.0f),
        vec2(-1.0f, 0.0f),
        vec2(0.0, -1.0f)
    );

    vec2 texel_size = 1.0f / textureSize(tex, 0);

    vec3 color = vec3(0.0f, 0.0f, 0.0f);
    float weight_sum;

    for (int i = 0; i < 4; i++) 
    {
		vec4 sampled_pos = SamplePositionAt(u_PositionTexture, txc + Kernel[i] * texel_size);

		if (sampled_pos.w <= 0.0f)
		{
			continue;
		}

        float sampled_depth = (sampled_pos.z); 
        float dweight = 1.0f / (abs(base_depth - sampled_depth) + 0.001f);

        float computed_weight = dweight;
        color.rgb += texture(tex, txc + Kernel[i] * texel_size).rgb * computed_weight;
        weight_sum += computed_weight;
    }

    color /= max(weight_sum, 0.2f);
    color = clamp(color, texture(tex, txc).rgb * 0.12f, vec3(1.0f));
    return color;
}

bool DetectAtEdge(in vec2 txc)
{
	vec2 TexelSize = 1.0f / textureSize(u_PositionTexture, 0);

	vec2 NeighbourhoodOffsets[8] = vec2[8]
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
		float T_at = SamplePositionAt(u_PositionTexture, txc + (NeighbourhoodOffsets[i] * TexelSize)).w;

		if (T_at <= 0.0f) 
		{
			return true;
		}
	}

	return false;
}

vec3 ToClipSpace(in vec3 pos)
{
	vec4 ViewSpace = u_ViewMatrix * vec4(pos, 1.0f);
    vec4 Projected = u_ProjectionMatrix * ViewSpace;
    Projected.xyz /= Projected.w;

    return Projected.xyz;
}

vec3 BilateralUpsample(sampler2D tex, vec2 txc, vec3 base_normal, float base_depth)
{
    const vec2 Kernel[4] = vec2[](
        vec2(0.0f, 1.0f),
        vec2(1.0f, 0.0f),
        vec2(-1.0f, 0.0f),
        vec2(0.0, -1.0f)
    );

    vec2 texel_size = 1.0f / textureSize(tex, 0);

    vec3 color = vec3(0.0f, 0.0f, 0.0f);
    float weight_sum;

    for (int i = 0; i < 4; i++) 
    {
        vec3 sampled_normal = SampleNormal(u_NormalTexture, txc + Kernel[i] * texel_size).xyz;
        float nweight = clamp(pow(abs(dot(sampled_normal, base_normal)), 8), 0.001f, 1.0f);
        float sampled_depth = texture(u_PositionTexture, txc + Kernel[i] * texel_size).x; 
        float dweight = 1.0f / (abs(base_depth - sampled_depth) + 0.01f);
		dweight = clamp(pow(dweight,32.0),0.01f,1.0f);
        float computed_weight = nweight * dweight;
        color.rgb += texture(tex, txc + Kernel[i] * texel_size).rgb * computed_weight;
        weight_sum += computed_weight;
    }

    color /= max(weight_sum, 0.01f);
    return color;
}

// Exponential Distance fog
vec3 ApplyFog(in vec3 InputColor, in float dist_to_point)
{
	float b = 0.00275f;
    float FogAmount = 1.0 - exp(-dist_to_point * b);
    vec3 FogColor  = vec3(0.5f, 0.6f, 0.7f);
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

vec3 BasicTonemap(vec3 color)
{
    float l = length(color);
    color = mix(color, color * 0.5f, l / (l + 1.0f));
    color = (color / sqrt(color * color + 1.0f));

    return color;
}



vec3 ColorSaturate(vec3 rgb, float adjustment)
{
    // Algorithm from Chapter 16 of OpenGL Shading Language
    const vec3 W = vec3(0.2125, 0.7154, 0.0721);
    vec3 intensity = vec3(dot(rgb, W));
    return mix(intensity, rgb, adjustment);
}

vec4 textureGood( sampler2D sam, vec2 uv )
{
    vec2 res = textureSize( sam, 0 );

    vec2 st = uv*res - 0.5;

    vec2 iuv = floor( st );
    vec2 fuv = fract( st );

    vec4 a = texture( sam, (iuv+vec2(0.5,0.5))/res );
    vec4 b = texture( sam, (iuv+vec2(1.5,0.5))/res );
    vec4 c = texture( sam, (iuv+vec2(0.5,1.5))/res );
    vec4 d = texture( sam, (iuv+vec2(1.5,1.5))/res );

    return mix( mix( a, b, fuv.x),
                mix( c, d, fuv.x), fuv.y );
}

float Brightness(vec3 c)
{
    return max(max(c.r, c.g), c.b);
}

vec3 HighQualityBloomUpsample(int lod) { 
	vec2 t = textureSize(u_BloomMips[lod], 0);
	float scale = 1.1f;
	vec4 d = t.xyxy * vec4(1, 1, -1, 0) * scale;
	vec2 uv = v_TexCoords;
    vec3 s = vec3(0.0f);
    s  = texture(u_BloomMips[lod], uv - d.xy).xyz;
    s += texture(u_BloomMips[lod], uv - d.wy).xyz * 2;
    s += texture(u_BloomMips[lod], uv - d.zy).xyz;
    s += texture(u_BloomMips[lod], uv + d.zw).xyz * 2;
    s += texture(u_BloomMips[lod], uv       ).xyz * 4;
    s += texture(u_BloomMips[lod], uv + d.xw).xyz * 2;
    s += texture(u_BloomMips[lod], uv + d.zy).xyz;
    s += texture(u_BloomMips[lod], uv + d.wy).xyz * 2;
    s += texture(u_BloomMips[lod], uv + d.xy).xyz;
    return s * (1.0 / 16);
}

//void SampleVolumetrics(float BaseLinearDepth, in vec3 BaseNormal, out vec3 TotalPointVL) 
//{
//	float TotalWeight = 0.0f;
//	vec2 TexelSize = 1.0f / textureSize(u_VolumetricsCompute, 0);
//	const float AtrousWeights[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);
//
//	TotalPointVL = vec3(0.0f);
//	vec3 GodRays = vec3(0.0f);
//
//	for (int x = -2 ; x <= 2 ; x++) {
//		for (int y = -2 ; y <= 2 ; y++) {
//			
//			vec2 SampleCoord = v_TexCoords + vec2(x,y) * TexelSize;
//			float LinearDepthAt = (texture(u_PositionTexture, SampleCoord).x);
//            if (LinearDepthAt <= 0.0f + 1e-2) { continue; }
//			vec3 NormalAt = SampleNormal(u_NormalTexture, SampleCoord.xy).xyz;
//			float DepthWeight = 1.0f / (abs(LinearDepthAt - BaseLinearDepth) + 0.01f);
//			DepthWeight = pow(DepthWeight, 12.0f);
//            DepthWeight = clamp(DepthWeight, 0.0f, 0.99f);
//			float NormalWeight = pow(abs(dot(NormalAt, BaseNormal)), 8.0f);
//			float Kernel = AtrousWeights[x + 2] * AtrousWeights[y + 2];
//			float Weight = Kernel * NormalWeight * DepthWeight;
//			Weight = max(Weight, 0.01f);
//			TotalPointVL += texture(u_VolumetricsCompute, SampleCoord).xyz * Weight;
//			TotalWeight += Weight;
//		}
//	}
//
//	TotalPointVL /= TotalWeight + 1e-4;
//}

vec4 textureSmooth(sampler2D t, vec2 x, vec2 textureSize)
{
    x *= vec2(textureSize);
    vec2 p = floor(x);
    vec2 f = fract(x);
    vec4 a = texture(t, (p                 ) / vec2(textureSize));
    vec4 b = texture(t, (p + vec2(1.0, 0.0)) / vec2(textureSize));
    vec4 c = texture(t, (p + vec2(0.0, 1.0)) / vec2(textureSize));
    vec4 d = texture(t, (p + vec2(1.0, 1.0)) / vec2(textureSize));
    return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

void FilmGrain(inout vec3 oc) 
{
	float Strength = u_FilmGrainStrength;
	vec3 NoiseColor = vec3(0.2001f, 0.804f, 1.02348f);
    vec3 Noise = vec3(hash2().x, hash2().y, hash2().x);
	oc += Noise * exp(-oc) * NoiseColor * 0.01f;
    oc *= mix(vec3(1.0f), Noise, Strength / 5.0f);
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
    
	// Barrel Distortion
    for (int i = 1; i <= 16; i++)
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

#define CONE_OVERLAP_SIMULATION 0.25

const mat3 HejlBurgessConeOverlapMatrix2Deg = mat3(
    mix(vec3(1.0, 0.0, 0.0), vec3(0.5595088340965042, 0.39845359892109633, 0.04203756698239944), vec3(CONE_OVERLAP_SIMULATION)),
    mix(vec3(0.0, 1.0, 0.0), vec3(0.43585871315661756, 0.5003841413971261, 0.06375714544625634), vec3(CONE_OVERLAP_SIMULATION)),
    mix(vec3(0.0, 0.0, 1.0), vec3(0.10997368482498855, 0.15247972169325025, 0.7375465934817612), vec3(CONE_OVERLAP_SIMULATION))
);

const mat3 InverseHejlBurgessConeOverlapMatrix2Deg = inverse(HejlBurgessConeOverlapMatrix2Deg);

vec3 TonemapHejlBurgess(in vec3 color) 
{
	color = (color * (6.2 * color + 0.5)) / (color * (6.2 * color + 1.7) + 0.06);
	return color;
}


vec3 vibrance(in vec3 color) 
{
	const vec3 lumacoeff_rec709 = vec3(0.2125, 0.7154, 0.0721);
    float lum = dot(color, lumacoeff_rec709);
    vec3 mask = (color - vec3(lum));
    mask = clamp(mask, 0.0, 1.0);
    float lum_mask = dot(lumacoeff_rec709, mask);
    lum_mask = 1.0 - lum_mask;
    return mix(vec3(lum), color, (1.0 + 0.2) * lum_mask);
}


vec4 SampleTextureCatmullRom(sampler2D tex, in vec2 uv); // catmull rom texture interp

 #define VOLUMETRIC_BICUBIC

vec4 texture_catmullrom(sampler2D tex, vec2 uv);
void main()
{
	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(u_Time) * 100.0f;



    float exposure = mix(u_LensFlare ? 3.77777f : 4.77777f, 1.25f, min(distance(-u_SunDirection.y, -1.0f), 0.99f));
	vec4 PositionAt = SamplePositionAt(u_PositionTexture, v_TexCoords).rgba;
	vec3 NormalAt = SampleNormal(u_NormalTexture, v_TexCoords).rgb;
	vec3 ClipSpaceAt = ToClipSpace(PositionAt.xyz);
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	float BayerDither = Bayer64(gl_FragCoord.xy);
	vec2 InverseVolRes = 1.0f / textureSize(u_VolumetricsCompute,0);
	vec3 PointVolumetrics = vec3(0.0f);
	
	if (u_PointVolumetricsToggled) {
	
		PointVolumetrics.xyz = texture_catmullrom(u_VolumetricsCompute, v_TexCoords+(Bayer128(gl_FragCoord.xy)*1.1*(1.0f/textureSize(u_VolumetricsCompute,0)))).xyz;
		//PointVolumetrics *= clamp(Bayer128(gl_FragCoord.xy),0.5f,1.0f);
		//PointVolumetrics.xyz = texture_catmullrom(u_VolumetricsCompute, v_TexCoords).xyz;
		//PointVolumetrics += BayerDither * (1.0f-exp(-GetLuminance(PointVolumetrics.xyz)));
	}



	if (PositionAt.w > 0.0f && !DetectAtEdge(v_TexCoords))
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
			float SampledSSAO = DepthOnlyBilateralUpsample(u_SSAOTexture, v_TexCoords, PositionAt.z).r;
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

			float RTAO = BilateralUpsample(u_RTAOTexture, v_TexCoords, NormalAt, PositionAt.w).r;
			RTAO = pow(RTAO, rtao_strength);
			RTAO = max(RTAO, 0.825f);

			InputColor = (RTAO * InputColor) + (0.01f * InputColor);
		}

		InputColor = ColorSaturate(InputColor, 1.1f);
		ColorGrading(InputColor);

		float fake_vol_multiplier = u_SunIsStronger ? 1.21f : 0.9f;

		if (u_GodRays)
		{
			float intensity = u_SunIsStronger ? 0.55f : 0.025f;
			float god_rays = DepthOnlyBilateralUpsample(u_VolumetricTexture, v_TexCoords, PositionAt.z).r * intensity * u_GodRaysStrength;
			vec3 ss_volumetric_color = u_SunIsStronger ? (vec3(250.0f, 250.0f, 230.0f) / 255.0f) : (vec3(96.0f, 192.0f, 255.0f) / 255.0f);
			InputColor += god_rays * ss_volumetric_color;
			fake_vol_multiplier = 0.3f;
		}

		if (u_SSGodRays)
		{
			float god_rays = GetScreenSpaceGodRays(PositionAt.xyz) * fake_vol_multiplier * u_GodRaysStrength;
			vec3 ss_volumetric_color = u_SunIsStronger ? (vec3(250.0f, 250.0f, 230.0f) / 255.0f) : (vec3(96.0f, 192.0f, 255.0f) / 255.0f);
			InputColor += god_rays * ss_volumetric_color;
		}

		InputColor += PointVolumetrics;

		if (u_ExponentialFog)
		{
			InputColor = ApplyFog(InputColor, PositionAt.w); 
			//InputColor = ApplyFog(InputColor, PositionAt.w, v_RayDirection, normalize(u_StrongerLightDirection));
		}

		o_Color = InputColor;
		float Exposure = mix(clamp(u_Exposure, 0.2f, 10.0f), clamp(u_Exposure - 1.25f, 0.2f, 10.0f), SunVisibility);
		if (u_PurkingeEffectStrength > 0.01f) {
			o_Color.xyz = PurkinjeEffect(o_Color.xyz);
		}

		if (u_HejlBurgess) {
			o_Color *= 0.275f * (Exposure * 0.5f);
			o_Color = TonemapHejlBurgess(o_Color * HejlBurgessConeOverlapMatrix2Deg) * InverseHejlBurgessConeOverlapMatrix2Deg;
			o_Color = vibrance(o_Color);
			o_Color = mix(vibrance(o_Color), o_Color, 0.5f);
		}

		else {
			o_Color = ACESFitted(vec4(o_Color, 1.0f), Exposure * 0.9f).rgb;
		}
	}

	else 
	{
		o_Color = u_ChromaticAberrationStrength <= 0.001f ? texture(u_FramebufferTexture, v_TexCoords).rgb : BasicChromaticAberation() ;
		o_Color += PointVolumetrics;
		o_Color *= clamp(clamp(u_Exposure * 0.5f, 0.0f, 10.0f) - 0.4256f, 0.0f, 10.0f);
		if (u_PurkingeEffectStrength>0.01f) {
			o_Color.xyz = PurkinjeEffect(o_Color.xyz);
		}
		o_Color = BasicTonemap(o_Color);
	}

	
	
	if (u_Bloom)
	{
		vec3 Bloom[4] = vec3[](vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f));

		// bicubic upsampling because the bloom textures are super low res
		// res : 0.25, 0.125, 0.0625, 0.03125
		// 16 total samples
		Bloom[0] += pow(textureBicubic(u_BloomMips[0], v_TexCoords).xyz, vec3(2.2f));
		Bloom[1] += pow(textureBicubic(u_BloomMips[1], v_TexCoords).xyz, vec3(2.2f)); 
		Bloom[2] += pow(textureBicubic(u_BloomMips[2], v_TexCoords).xyz, vec3(2.2f)); 
		Bloom[3] += pow(textureBicubic(u_BloomMips[3], v_TexCoords).xyz, vec3(2.2f)); 

		vec3 TotalBloom = vec3(0.0f);
		float AmplificationFactor = 1.350f;
		
		// fit the bloom to a curve which amplifies it based on some factor --->

		float Weights[4] = float[4](0.5f, 0.25f, 0.25f, 0.25f);
		TotalBloom = (pow(Bloom[0] ,vec3(1.0f / AmplificationFactor)) * Weights[0]) + TotalBloom;
		TotalBloom = (pow(Bloom[1] ,vec3(1.0f / AmplificationFactor)) * Weights[1])  + TotalBloom;
		TotalBloom = (pow(Bloom[2] ,vec3(1.0f / AmplificationFactor)) * Weights[2])  + TotalBloom;
		TotalBloom = (pow(Bloom[3] ,vec3(1.0f / AmplificationFactor)) * Weights[3]) + TotalBloom;
		o_Color += TotalBloom;
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

	if (u_FilmGrain) {
		FilmGrain(o_Color);
	}
	


	// used to test out ways to animate blue noise

	///ivec2 TxS = ivec2(textureSize(u_BlueNoise, 0).xy);
	///vec3 Hash = texelFetch(u_BlueNoise, ivec2(gl_FragCoord.xy*2.0)%TxS, 0).xyz;//texture(u_BlueNoise, v_TexCoords * (u_Dimensions / vec2(TxS / 2.0f))).xy;
	///const float GoldenRatio = 1.61803398875;
	///vec3 Xi = mod(Hash + GoldenRatio * (u_CurrentFrame % 240), 1.0f);
	///o_Color = vec3(Xi);

	///vec3 Hash;
	///int n = u_CurrentFrame%1024;
	///vec2 off = fract(vec2(n*12664745, n*9560333)/16777216.0) * 1024.0;
	///ivec2 TextureSize = textureSize(u_BlueNoise, 0);
	///ivec2 SampleTexelLoc = ivec2(gl_FragCoord.xy + ivec2(floor(off))) % TextureSize;
	///Hash = texelFetch(u_BlueNoise, SampleTexelLoc, 0).xyz;
	///vec3 Xi = Hash;
	///o_Color = Xi;

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

vec4 sampleLevel0(sampler2D tex, vec2 uv )
{
    return texture( tex, uv, -10.0 );
}

// note: entirely stolen from https://gist.github.com/TheRealMJP/c83b8c0f46b63f3a88a5986f4fa982b1
//
// Samples a texture with Catmull-Rom filtering, using 9 texture fetches instead of 16.
// See http://vec3.ca/bicubic-filtering-in-fewer-taps/ for more details
vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv)
{
	vec2 texSize = textureSize(tex, 0);

    // We're going to sample a a 4x4 grid of texels surrounding the target UV coordinate. We'll do this by rounding
    // down the sample location to get the exact center of our "starting" texel. The starting texel will be at
    // location [1, 1] in the grid, where [0, 0] is the top left corner.
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5) + 0.5;

    // Compute the fractional offset from our starting texel to our original sample location, which we'll
    // feed into the Catmull-Rom spline function to get our filter weights.
    vec2 f = samplePos - texPos1;

    // Compute the Catmull-Rom weights using the fractional offset that we calculated earlier.
    // These equations are pre-expanded based on our knowledge of where the texels will be located,
    // which lets us avoid having to evaluate a piece-wise function.
    vec2 w0 = f * ( -0.5 + f * (1.0 - 0.5*f));
    vec2 w1 = 1.0 + f * f * (-2.5 + 1.5*f);
    vec2 w2 = f * ( 0.5 + f * (2.0 - 1.5*f) );
    vec2 w3 = f * f * (-0.5 + 0.5 * f);
    
    // Work out weighting factors and sampling offsets that will let us use bilinear filtering to
    // simultaneously evaluate the middle 2 samples from the 4x4 grid.
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / w12;

    // Compute the final UV coordinates we'll use for sampling the texture
    vec2 texPos0 = texPos1 - vec2(1.0);
    vec2 texPos3 = texPos1 + vec2(2.0);
    vec2 texPos12 = texPos1 + offset12;

    texPos0 /= texSize;
    texPos3 /= texSize;
    texPos12 /= texSize;

    vec4 result = vec4(0.0);
    result += sampleLevel0(tex, vec2(texPos0.x,  texPos0.y)) * w0.x * w0.y;
    result += sampleLevel0(tex, vec2(texPos12.x, texPos0.y)) * w12.x * w0.y;
    result += sampleLevel0(tex, vec2(texPos3.x,  texPos0.y)) * w3.x * w0.y;

    result += sampleLevel0(tex, vec2(texPos0.x,  texPos12.y)) * w0.x * w12.y;
    result += sampleLevel0(tex, vec2(texPos12.x, texPos12.y)) * w12.x * w12.y;
    result += sampleLevel0(tex, vec2(texPos3.x,  texPos12.y)) * w3.x * w12.y;

    result += sampleLevel0(tex, vec2(texPos0.x,  texPos3.y)) * w0.x * w3.y;
    result += sampleLevel0(tex, vec2(texPos12.x, texPos3.y)) * w12.x * w3.y;
    result += sampleLevel0(tex, vec2(texPos3.x,  texPos3.y)) * w3.x * w3.y;

    return result;
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