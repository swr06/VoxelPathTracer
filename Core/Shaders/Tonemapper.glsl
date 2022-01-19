#version 330 core
#define CONE_OVERLAP_SIMULATION 0.25

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_DOFTexture;
uniform sampler2D u_Depth;
uniform bool u_FilmGrain;
uniform bool u_HejlBurgess;
uniform bool u_DOF;
uniform float u_FilmGrainStrength;
uniform float u_Time;
uniform float u_Exposure;
uniform vec2 u_CameraPlanes;
uniform float u_FocalDepthTemporal;
uniform float u_COCScale;

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

// Hejl Burgess ->
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

// ACES ->
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
    Color.rgb *= Exposure;
    Color.rgb = ACESInputMat * Color.rgb;
    Color.rgb = RRTAndODTFit(Color.rgb);
    Color.rgb = ACESOutputMat * Color.rgb;
    return Color;
}

vec2 g_TexCoords;

// Basic cubic distortion 
vec2 CubicDistort(vec2 uv) {
	float k = -1.0f * 0.085f, kcube = 0.5f * 0.4f;
    float r2 = (uv.x - 0.5f) * (uv.x - 0.5f) + (uv.y - 0.5f) * (uv.y - 0.5f);
    float f = 0.0f;
    if (kcube == 0.0f) f = 1.0f + r2 * k;
    else f = 1.0f + r2 * (k + kcube * sqrt(r2));
    float x = f * (uv.x - 0.5f) + 0.5f;
    float y = f * (uv.y - 0.5f) + 0.5f;
	return vec2(x, y);
}

vec4 textureBicubic(sampler2D sampler, vec2 texCoords);

// https://en.wikipedia.org/wiki/Circle_of_confusion
float GetCircleOfConfusion(float Depth, float CenterDepth, float Scale) 
{
	float CircleOfConfusion = abs(Depth - CenterDepth) / 0.6f;
	return CircleOfConfusion / (1.0f / Scale + CircleOfConfusion);
}

float DelinearizeDepth(float d) 
{
	float near = u_CameraPlanes.x;
	float far = u_CameraPlanes.y;
    d = d < 0.0f ? d : d + near;
	d = d < 0.0f ? far - 0.01f : d;
	float Delinearized = -((near + far) * d - (2.0f * near)) / ((near - far) * d);
	return Delinearized;
}


// Main ->
void main() {

	HASH2SEED = (g_TexCoords.x * g_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(u_Time) * 100.0f;

	// Distort ->
	g_TexCoords = v_TexCoords;
	//g_TexCoords = CubicDistort(g_TexCoords);


    // Linear depth ->
    float BaseDepth = texture(u_Depth, g_TexCoords).x;


	vec3 BaseColor = texture(u_InputTexture, g_TexCoords).xyz;
    o_Color = BaseColor;


    // Mix DOF 
    if (u_DOF) {
        
        // CoC is calculated using delinearized depth
        float DelinearizedCenterDepth = DelinearizeDepth(BaseDepth);
        
        vec3 DOFFetch = textureBicubic(u_DOFTexture, g_TexCoords).xyz;

	    float CoC = GetCircleOfConfusion(DelinearizedCenterDepth, u_FocalDepthTemporal, u_COCScale);
        float MixFactor = CoC * 2050.0f;

        MixFactor = clamp(MixFactor, 0.0f, 1.0f);

        o_Color = mix(BaseColor, DOFFetch, MixFactor);
    }


	float Exposure = u_Exposure * 0.9f * 0.58f;

	// Tonemap ->
	const bool TONE_MAP = true;
	if (TONE_MAP) {


		if (u_HejlBurgess) {
			o_Color *= 0.275f * (Exposure * 0.5f);
			o_Color = TonemapHejlBurgess(o_Color * HejlBurgessConeOverlapMatrix2Deg) * InverseHejlBurgessConeOverlapMatrix2Deg;
			o_Color = vibrance(o_Color);
			o_Color = mix(vibrance(o_Color), o_Color, 0.5f);
		}

		else {
			o_Color = ACESFitted(vec4(o_Color, 1.0f), Exposure).rgb;
			o_Color = mix(vibrance(o_Color), o_Color, 0.6f);
		}
	}

	if (u_FilmGrain) {
		FilmGrain(o_Color);
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