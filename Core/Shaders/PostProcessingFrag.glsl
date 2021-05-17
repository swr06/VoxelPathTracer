#version 330 core

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

uniform vec3 u_SunDirection;
uniform vec3 u_StrongerLightDirection;
uniform vec2 u_Dimensions;

uniform bool u_SunIsStronger;

uniform bool u_LensFlare = true;
uniform bool u_GodRays = true;
uniform bool u_SSAO = false;

uniform sampler2D u_FramebufferTexture;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;

uniform sampler2D u_BlueNoise;
uniform sampler2D u_SSAOTexture;

uniform mat4 u_ProjectionMatrix;
uniform mat4 u_ViewMatrix;

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
    const int SAMPLES = 8;
	float dither = texture(u_BlueNoise, v_TexCoords * (u_Dimensions / vec2(256.0f))).r;

    for (int i = 0; i < SAMPLES; i++)
    {
        float scale = (1.0f - (float(i) / float(SAMPLES))) + dither / float(SAMPLES);

        vec2 coord = (v_TexCoords - SunScreenSpacePosition) * scale + SunScreenSpacePosition;
        coord = clamp(coord, 0.001f, 0.999f);

        float is_sky_at = texture(u_PositionTexture, coord).w > 0.0f ? 0.0f : 1.0f;

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
		vec4 sampled_pos = texture(u_PositionTexture, txc + Kernel[i] * texel_size);

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
		float T_at = texture(u_PositionTexture, txc + (NeighbourhoodOffsets[i] * TexelSize)).w;

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

void main()
{
    float exposure = mix(u_LensFlare ? 3.77777f : 4.77777f, 1.25f, min(distance(-u_SunDirection.y, -1.0f), 0.99f));
	vec4 PositionAt = texture(u_PositionTexture, v_TexCoords).rgba;
	bool AtEdge = DetectAtEdge(v_TexCoords);
	vec3 ClipSpaceAt = ToClipSpace(PositionAt.xyz);

	if (PositionAt.w > 0.0f && (!AtEdge))
	{
		vec3 InputColor;
		vec3 Sharpened = sharpen(u_FramebufferTexture, v_TexCoords);
		InputColor = texture(u_FramebufferTexture, v_TexCoords).rgb;

		if (u_SSAO && ((1.0f - ClipSpaceAt.z) > 0.006f))
		{
			const float ssao_strength = 4.0f;
			float SampledSSAO = DepthOnlyBilateralUpsample(u_SSAOTexture, v_TexCoords, PositionAt.z).r;
			//float SampledSSAO = texture(u_SSAOTexture, v_TexCoords).r;
			float SSAO = pow(SampledSSAO, ssao_strength);
			SSAO = clamp(SSAO, 0.0001, 1.0f);
			InputColor *= ssao_strength - 1.1f;
			InputColor *= SSAO;
		}

		ColorGrading(InputColor);
		ColorSaturation(InputColor);

		if (u_GodRays)
		{
			float god_rays = GetScreenSpaceGodRays(PositionAt.xyz);
			vec3 ss_volumetric_color = u_SunIsStronger ? (vec3(189.0f, 200.0f, 129.0f) / 255.0f) : (vec3(96.0f, 192.0f, 255.0f) / 255.0f);
			InputColor += ss_volumetric_color * god_rays;
		}

		InputColor = ACESFitted(vec4(InputColor, 1.0f), exposure + 0.01f).rgb;
		o_Color = InputColor;
	}

	else 
	{
		o_Color = texture(u_FramebufferTexture, v_TexCoords).rgb;
	}

	if (u_LensFlare && u_SunIsStronger)
	{
		vec2 SunScreenSpacePosition = WorldToScreen(u_SunDirection * 10000.0f) - 0.5f; 
		SunScreenSpacePosition.x *= u_Dimensions.x / u_Dimensions.y;
		vec2 LensFlareCoord = v_TexCoords - 0.5f;
		LensFlareCoord.x *= u_Dimensions.x / u_Dimensions.y;
		
		vec3 LensFlare = vec3(1.6f, 1.2f, 1.0f) * lensflare(LensFlareCoord, SunScreenSpacePosition);
		LensFlare = clamp(LensFlare, 0.02f, 0.999f);
		o_Color += LensFlare;
	}

}