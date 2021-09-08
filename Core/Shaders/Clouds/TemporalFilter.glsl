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

    for(int x = -3; x <= 3; x++) 
	{
        for(int y = -3; y <= 3; y++) 
		{
            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize).rgb; 
            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

    return clipAABB(Color, MinColor, MaxColor);
}



vec4 textureBicubic(sampler2D sampler, vec2 texCoords);

void main()
{
	vec4 CurrentColor = textureBicubic(u_CurrentColorTexture, v_TexCoords).rgba;

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
		vec4 PrevColor = textureBicubic(u_PreviousColorTexture, PreviousCoord).rgba;

		if (u_Clamp) {
			PrevColor.xyz = ClampColor(PrevColor.xyz);
		}
		float T_Base = texture(u_CurrentPositionData, v_TexCoords).r;
		float T_Prev = texture(u_PrevPositionData, PreviousCoord).r;

		bool OcclusionValidity = (T_Base > 0.0f) == (T_Prev > 0.0f);

		float ReprojectBias = 0.006458f;
		if(PreviousCoord.x > 0.0 + ReprojectBias && PreviousCoord.x < 1.0 - ReprojectBias
			&& PreviousCoord.y > 0.0 + ReprojectBias && PreviousCoord.y < 1.0 - ReprojectBias && 
			OcclusionValidity)
		{
			o_Color = mix(CurrentColor.xyzw, PrevColor.xyzw, 0.925f);
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