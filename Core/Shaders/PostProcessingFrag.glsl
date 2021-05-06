#version 330 core

layout(location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;

uniform sampler2D u_FramebufferTexture;
uniform sampler2D u_PositionTexture;

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
	sum += -1. * smoothfilter(tex, coords + vec2( -1.0 * dx , 0.0 * dy));
	sum += -1. * smoothfilter(tex, coords + vec2( 0.0 * dx , -1.0 * dy));
	sum += 5. * smoothfilter(tex, coords + vec2( 0.0 * dx , 0.0 * dy));
	sum += -1. * smoothfilter(tex, coords + vec2( 0.0 * dx , 1.0 * dy));
	sum += -1. * smoothfilter(tex, coords + vec2( 1.0 * dx , 0.0 * dy));
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

void main()
{
	vec2 TexSize = textureSize(u_PositionTexture, 0);
    float PixelDepth1 = texture(u_PositionTexture, clamp(v_TexCoords + vec2(0.0f, 1.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;
    float PixelDepth2 = texture(u_PositionTexture, clamp(v_TexCoords + vec2(0.0f, -1.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;
    float PixelDepth3 = texture(u_PositionTexture, clamp(v_TexCoords + vec2(1.0f, 0.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;
    float PixelDepth4 = texture(u_PositionTexture, clamp(v_TexCoords + vec2(-1.0f, 0.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;

	if (texture(u_PositionTexture, v_TexCoords).w > 0.0f && PixelDepth1 > 0.0f && PixelDepth2 > 0.0f && PixelDepth3 > 0.0f && PixelDepth4 > 0.0f)
	{
		vec3 SharpenedColor = sharpen(u_FramebufferTexture, v_TexCoords);
		SharpenedColor = ACESFitted(vec4(SharpenedColor, 1.0f), 4.5f).rgb;
		Vignette(SharpenedColor);
		o_Color = SharpenedColor;
	}

	else 
	{
		o_Color = texture(u_FramebufferTexture, v_TexCoords).rgb;
	}
}