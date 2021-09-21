#version 330 core
#define DENOISE
//#define BE_USELESS

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_CurrentPositionData;
uniform sampler2D u_PrevPositionData;

uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_MixModifier = 0.8;
uniform float u_Time;

uniform bool u_Clamp;

uniform vec3 u_CurrentPosition;
uniform vec3 u_PreviousPosition;

vec2 Reprojection(vec3 pos) 
{
	vec3 WorldPos = pos;
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xy;
}

vec2 ProjectCurrent(vec3 pos) 
{
	vec3 WorldPos = pos;
	vec4 ProjectedPosition = u_Projection * u_View * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xy;
}

float FastDist(vec3 p1, vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

float GetLuminance(vec3 color)
{
	return dot(color, vec3(0.299f, 0.587f, 0.114f));
}

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 GetBlurredColor() {

	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);

	vec4 TotalColor = vec4(0.0f);
	float Samples = 0.0f;

	// 6x6 taps
	for (int x = -3 ; x <= 3; x++) {
		for (int y = -3 ; y <= 3 ; y++) {
			vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;
			TotalColor += texture(u_CurrentColorTexture, SampleCoord);
			Samples += 1.0f;
		}
	}

	TotalColor /= Samples;
	return TotalColor;
}

vec3 clipAABB(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

vec3 ClampColor(vec3 Color) 
{
    vec3 MinColor = vec3(100.0);
	vec3 MaxColor = vec3(-100.0); 
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture,0);

    for(int x = -1; x <= 1; x++) 
	{
        for(int y = -2; y <= 2; y++) 
		{
            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize).rgb; 
            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

    float Bias = 0.08f;
    return clamp(Color, MinColor-Bias , MaxColor+Bias);
}



vec4 textureBicubic(sampler2D sampler, vec2 texCoords);
vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv, vec2 texSize );
vec4 texture_catmullrom(sampler2D tex, vec2 uv);

void main()
{
	vec4 CurrentColor = texture_catmullrom(u_CurrentColorTexture, v_TexCoords).rgba;
	//vec4 CurrentColor = texture(u_CurrentColorTexture, v_TexCoords).rgba;
	//vec4 CurrentColor = SampleTextureCatmullRom(u_CurrentColorTexture, v_TexCoords, textureSize(u_CurrentColorTexture,0)).rgba;


	#ifdef BE_USELESS
	o_Color = CurrentColor;
	return;
	#endif

	vec3 PlayerPosition = u_InverseView[3].xyz;

	if (CurrentColor.w > -0.5f)
	{
		// Reproject clouds by assuming it lies in a plane some x units in front of the camera :
		vec3 CurrentVirtualPosition = PlayerPosition + normalize(GetRayDirectionAt(v_TexCoords)) * 800.0f;
		vec2 ProjectedCurrent = ProjectCurrent(CurrentVirtualPosition);
		vec2 ProjectedPrevious = Reprojection(CurrentVirtualPosition.xyz);
		
		vec2 PreviousCoord = ProjectedPrevious * 0.5f + 0.5f;
		ProjectedCurrent = ProjectedCurrent * 0.5f + 0.5f;
		//vec4 PrevColor = SampleTextureCatmullRom(u_PreviousColorTexture, PreviousCoord, textureSize(u_PreviousColorTexture,0)).rgba;
		vec4 PrevColor = texture_catmullrom(u_PreviousColorTexture, PreviousCoord).rgba;

		if (u_Clamp) {
			PrevColor.xyz = ClampColor(PrevColor.xyz);
		}

		float T_Base = texture(u_CurrentPositionData, v_TexCoords).r;
		float T_Prev = texture(u_PrevPositionData, PreviousCoord).r;
		const float Thresh = sqrt(2.0f) * 6.942069420f;
		bool OcclusionValidity = (T_Base > 0.0f) == (T_Prev > 0.0f);

		float ReprojectBias = 0.011f;
		if(PreviousCoord.x > 0.0 + ReprojectBias && PreviousCoord.x < 1.0 - ReprojectBias
			&& PreviousCoord.y > 0.0 + ReprojectBias && PreviousCoord.y < 1.0 - ReprojectBias && 
			OcclusionValidity)
		{
			o_Color = mix(CurrentColor.xyzw, PrevColor.xyzw, 0.94750f);
		}

		else 
		{
		#ifdef DENOISE
			o_Color = GetBlurredColor();
		#else 
			o_Color = CurrentColor;
		#endif
		}
	}

	else 
	{
		#ifdef DENOISE 
			o_Color = GetBlurredColor();
		#else 
			o_Color = CurrentColor;
		#endif
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

float sqr(float x) { return x*x; }
vec2 sqr(vec2 x) { return x*x; }

vec4 sampleLevel0(sampler2D tex, vec2 uv )
{
    return texture( tex, uv, -10.0 );
}

vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv, vec2 texSize )
{
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5) + 0.5;
    vec2 f = samplePos - texPos1;
    vec2 w0 = f * ( -0.5 + f * (1.0 - 0.5*f));
    vec2 w1 = 1.0 + f * f * (-2.5 + 1.5*f);
    vec2 w2 = f * ( 0.5 + f * (2.0 - 1.5*f) );
    vec2 w3 = f * f * (-0.5 + 0.5 * f);
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / (w1 + w2);
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
        col    += textureLod(tex, vec2(uv0.x, uv0.y), 0)*w0.x*w0.y;
        col    += textureLod(tex, vec2(uv12.x, uv0.y), 0)*w12.x*w0.y;
        col    += textureLod(tex, vec2(uv3.x, uv0.y), 0)*w3.x*w0.y;

        col    += textureLod(tex, vec2(uv0.x, uv12.y), 0)*w0.x*w12.y;
        col    += textureLod(tex, vec2(uv12.x, uv12.y), 0)*w12.x*w12.y;
        col    += textureLod(tex, vec2(uv3.x, uv12.y), 0)*w3.x*w12.y;

        col    += textureLod(tex, vec2(uv0.x, uv3.y), 0)*w0.x*w3.y;
        col    += textureLod(tex, vec2(uv12.x, uv3.y), 0)*w12.x*w3.y;
        col    += textureLod(tex, vec2(uv3.x, uv3.y), 0)*w3.x*w3.y;

    return clamp(col, 0.0, 65535.0);
}