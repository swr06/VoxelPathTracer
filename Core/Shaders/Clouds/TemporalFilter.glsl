#version 430 core
//#define DENOISE
//#define BE_USELESS

layout (location = 0) out vec4 o_Color;

// Equi rectangular projected clouds ->
layout(rgba8, binding = 4) uniform image2D o_EquiangularProjection;


// Bayer dither 
float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))


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
uniform bool u_SpatialUpscale;

uniform vec3 u_CurrentPosition;
uniform vec3 u_PreviousPosition;

uniform bool u_UpdateProjection;

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

vec4 SpatialBasic() {

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
        for(int y = -1; y <= 1; y++) 
		{
            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize * 1.1f).rgb; 
            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

    const float Bias = 0.05f;
    return clamp(Color, MinColor-Bias , MaxColor+Bias);
}

float Luminance(vec3 x)
{
	return dot(x, vec3(0.2125f, 0.7154f, 0.0721f));
}

vec4 texture_catmullrom(sampler2D tex, vec2 uv);

vec4 SpatialUpscaleCloud(sampler2D samp, vec2 txc) 
{
    vec2 TexelSize = 1.0f/textureSize(samp,0).xy;
    vec2 Kernel[4] = vec2[4](vec2(1.0f, 0.0f), vec2(0.0f, 1.0f), vec2(-1.0f, 0.0f), vec2(0.0f, -1.0f));
    vec4 BaseData = texture(samp, txc);
    float BaseLuminance = Luminance(BaseData.xyz);
    vec4 Blurred = vec4(0.0f);
    float Descrepancy = 0.0f;
    int InvalidSamples = 0;

    for (int x = 0 ; x < 4 ; x++) {
        vec2 Coord = txc+Kernel[x]*TexelSize;
        vec4 Sample = texture(samp, Coord).xyzw;
        float SampleLuma = Luminance(Sample.xyz);
        Blurred += Sample;
        float E = abs(BaseLuminance - SampleLuma);
        Descrepancy += E;
        
        if (E > 0.1f) 
        {
            InvalidSamples++;
        }
    }

    Blurred /= 4.0f;
    
    bool AntiChecker = InvalidSamples >= 3;
    return AntiChecker ? Blurred : texture_catmullrom(samp, txc).xyzw;
}

float pow5(float x) {
    return x*x*x*x*x;
}

float pow8(float x) {
    return x*x*x*x*x*x*x*x;
}

float pow12(float x) {
    return x*x*x*x*x*x*x*x*x*x*x*x;
}

vec4 SpatialFilterCloud(vec4 CurrentColor) {
    
    float BaseLuminance = Luminance(CurrentColor.xyz);
    
    vec4 TotalColor = vec4(0.0f);
    vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0).xy;
    const float Scale = 1.05f;
    float TotalWeight = 0.0f;

    // 36 taps.
    for (int x = -3 ; x <= 3 ; x++)
    {
        for (int y = -3 ; y <= 3 ; y++) 
        {
            vec2 SampleTxc = v_TexCoords+vec2(x,y)*Scale*TexelSize;
            //float XWeight = AtrousWeights[abs(x)];
	        //float YWeight = AtrousWeights[abs(y)];
            //float Kernel = clamp(XWeight*YWeight, 0.0f, 1.0f);
            vec4 Sample = texture(u_CurrentColorTexture, SampleTxc);
            float LumaAt = Luminance(Sample.xyz);
            float LuminanceError = clamp(1.0f - clamp(abs(LumaAt - BaseLuminance) / 3.0f, 0.0f, 1.0f), 0.0f, 1.0f);
            float Weight = clamp(pow12(LuminanceError), 0.0f, 1.0f); //Kernel * clamp(pow(LuminanceError, 3.5f), 0.0f, 1.0f);
            Weight = clamp(Weight, 0.01f, 1.0f);
            TotalColor += Sample*Weight;
            TotalWeight += Weight;
        }
    }

    return TotalColor / max(TotalWeight, 0.01f);
}


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

// Luminosity weights ->
// Effective as fuck.

vec4 SpatialUpscaleCloud2(sampler2D samp, vec2 txc, bool LargeKernel, bool StrictAssExponent, out vec4 MinData, out vec4 MaxData) 
{
    vec2 TexelSize = 1.0f/textureSize(samp,0).xy;
    vec4 BaseData = texture(samp, txc);
    float BaseLuminance = Luminance(BaseData.xyz);
    float TotalWeight = 1.0f;
    vec4 TotalData = BaseData;

    int KernelSize = LargeKernel ? 2 : 1;

    float LumaExp = StrictAssExponent ? 1200.0f : 55.0f;

    MinData = vec4(1000.0f);
    MaxData = vec4(-1000.0f);

    for (int x = -KernelSize ; x <= KernelSize ; x++) 
    {
        for (int y = -KernelSize ; y <= KernelSize ; y++) {

            if (x == 0 && y == 0) { continue ; }
            vec2 SampleCoord = txc+vec2(x,y)*TexelSize;
            
            if (SampleCoord != clamp(SampleCoord, 0.001f, 0.999f)) {
                continue;
            }

            vec4 SampleData = BetterTexture(samp, SampleCoord).xyzw;
            float LumaAt = Luminance(SampleData.xyz);
            float LumaError = 1.0f - clamp(abs(LumaAt - BaseLuminance) / 4.0f, 0.0f, 1.0f);
            float Weight = pow(abs(LumaError), LumaExp);
            TotalData += SampleData * Weight;
            TotalWeight += Weight;

            MinData = min(SampleData, MinData);
            MaxData = max(SampleData, MaxData);
        }
    }

    MinData = min(BaseData, MinData);
    MaxData = max(BaseData, MaxData);

    TotalData /= TotalWeight;
    return TotalData;
}




vec4 textureBicubic(sampler2D sampler, vec2 texCoords);
vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv, vec2 texSize );
vec4 texture2D_bicubic(sampler2D tex, vec2 uv);

bool SampleValid(in vec2 txc)
{
    const ivec2 Kernel[8] = ivec2[8]
	(
		ivec2(-1.0, -1.0),
		ivec2( 0.0, -1.0),
		ivec2( 1.0, -1.0),
		ivec2(-1.0,  0.0),
		ivec2( 1.0,  0.0),
		ivec2(-1.0,  1.0),
		ivec2( 0.0,  1.0),
		ivec2( 1.0,  1.0)
	);

	ivec2 basecoord = ivec2(floor(txc*textureSize(u_CurrentPositionData,0).xy));

    for (int i = 0 ; i < 8 ; i++)
    {
		ivec2 samplecoord = basecoord+Kernel[i];
        //if (texture(u_PositionTex, txc + Kernel[i] * TexelSize * 1.1f).r <= 0.0f)
		if (texelFetch(u_CurrentPositionData, samplecoord,0).x <= 0.0f)
        {
            return true;
        }
    }

    return false;
}

const float PI = 3.141592653f;


vec2 ProjectDirection(vec3 Direction) 
{
    vec2 TextureSize = vec2(imageSize(o_EquiangularProjection).xy);
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

 

vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv);

void main()
{
	vec4 CurrentColor = vec4(0.0f);

    vec4 MinData = vec4(-1000.0f), MaxData = vec4(1000.0f);

    if (u_SpatialUpscale) {
        CurrentColor = SpatialUpscaleCloud2(u_CurrentColorTexture, v_TexCoords, false, false, MinData, MaxData).xyzw;
    }

    else 
    {
        CurrentColor = texture_catmullrom(u_CurrentColorTexture, v_TexCoords).xyzw;
    }


	#ifdef BE_USELESS
	o_Color = CurrentColor;
	return;
	#endif

	vec3 PlayerPosition = u_InverseView[3].xyz;

    vec3 RayDirection = normalize(GetRayDirectionAt(v_TexCoords));

    bool SampleIsValid = true;//SampleValid(v_TexCoords);

	if (CurrentColor.w > -0.5f && true)
	{
		// Reproject clouds by assuming it lies in a plane some x units in front of the camera :

		vec3 CurrentVirtualPosition = PlayerPosition + RayDirection * 500.0f;
		vec2 ProjectedCurrent = ProjectCurrent(CurrentVirtualPosition);
		vec2 ProjectedPrevious = Reprojection(CurrentVirtualPosition.xyz);
		
		vec2 PreviousCoord = ProjectedPrevious * 0.5f + 0.5f;
		ProjectedCurrent = ProjectedCurrent * 0.5f + 0.5f;
        vec4 PrevColor;

        
        if (u_SpatialUpscale) 
        {
            float bayerdither = bayer256(gl_FragCoord.xy) / 2.8f;
            PreviousCoord += vec2(bayerdither, bayerdither * 0.5f) * (1.0f / textureSize(u_PreviousColorTexture, 0).xy);
            vec4 t1,t2;
            PrevColor = SampleTextureCatmullRom(u_PreviousColorTexture, PreviousCoord).xyzw;
            //PrevColor = FastCatmullRom(u_PreviousColorTexture, PreviousCoord, 0.5).xyzw;

            //PrevColor = SpatialUpscaleCloud2(u_PreviousColorTexture, PreviousCoord, true, true,t1,t2).xyzw;
           // PrevColor = textureBicubic(u_PreviousColorTexture, PreviousCoord).xyzw;
          //  PrevColor = texelFetch(u_PreviousColorTexture, ivec2(PreviousCoord*textureSize(u_PreviousColorTexture,0)), 0).xyzw;
        }

        else {
            PrevColor = SampleTextureCatmullRom(u_PreviousColorTexture, PreviousCoord).xyzw;
        }

		if (u_Clamp) {
			//PrevColor.xyz = ClampColor(PrevColor.xyz);
            PrevColor = clamp(PrevColor, MinData - 0.35f, MaxData + 0.35f);
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
			o_Color = mix(CurrentColor.xyzw, PrevColor.xyzw, 0.925f); 
		}

		else 
		{
		#ifdef DENOISE
        if (!u_SmartUpscale) {
			o_Color = SpatialBasic();
        }

        else {
            o_Color = SpatialFilterCloud(CurrentColor);
        }
		#else 
			o_Color = CurrentColor;
		#endif
		}
	}

	else 
	{

		#ifdef DENOISE 
	    if (!u_SmartUpscale) {
			o_Color = SpatialBasic();
        }

        else {
            o_Color = SpatialFilterCloud(CurrentColor);
        }

		#else 
			o_Color = CurrentColor;
		#endif
	}

    // Project and store in equiangular texture for gi ->

    const bool StoreEquiangular = true;

    if (StoreEquiangular && u_UpdateProjection) {

        vec2 UVProjection = ProjectDirection(RayDirection);

        if (UVProjection == clamp(UVProjection, 0.001f, 0.999f)) {

            ivec2 PixelProjection;
            PixelProjection = ivec2(floor(UVProjection * vec2(imageSize(o_EquiangularProjection).xy)));
            imageStore(o_EquiangularProjection, PixelProjection, vec4(o_Color.xyz, clamp(o_Color.w, 0.5f, 1.0f)));
        }
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
    return texture( tex, uv );
}

vec4 SampleTextureCatmullRom(sampler2D tex, vec2 uv)
{
    vec2 texSize = textureSize(tex, 0);
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



// Somewhat better bicubic interpolation 
// used for testing

float w0(float a)
{
    return (1.0/6.0)*(a*(a*(-a + 3.0) - 3.0) + 1.0);
}

float w1(float a)
{
    return (1.0/6.0)*(a*a*(3.0*a - 6.0) + 4.0);
}

float w2(float a)
{
    return (1.0/6.0)*(a*(a*(-3.0*a + 3.0) + 3.0) + 1.0);
}

float w3(float a)
{
    return (1.0/6.0)*(a*a*a);
}

float g0(float a)
{
    return w0(a) + w1(a);
}

float g1(float a)
{
    return w2(a) + w3(a);
}

float h0(float a)
{
    return -1.0 + w1(a) / (w0(a) + w1(a));
}

float h1(float a)
{
    return 1.0 + w3(a) / (w2(a) + w3(a));
}

vec4 texture2D_bicubic(sampler2D tex, vec2 uv)
{
    vec2 texel_size = 1.0f/textureSize(tex,0);
	vec4 texelSize = vec4(texel_size,1.0/texel_size);
	uv = uv*texelSize.zw;
	vec2 iuv = floor( uv );
	vec2 fuv = fract( uv );

    float g0x = g0(fuv.x);
    float g1x = g1(fuv.x);
    float h0x = h0(fuv.x);
    float h1x = h1(fuv.x);
    float h0y = h0(fuv.y);
    float h1y = h1(fuv.y);

	vec2 p0 = (vec2(iuv.x + h0x, iuv.y + h0y) - 0.5) * texelSize.xy;
	vec2 p1 = (vec2(iuv.x + h1x, iuv.y + h0y) - 0.5) * texelSize.xy;
	vec2 p2 = (vec2(iuv.x + h0x, iuv.y + h1y) - 0.5) * texelSize.xy;
	vec2 p3 = (vec2(iuv.x + h1x, iuv.y + h1y) - 0.5) * texelSize.xy;

    return g0(fuv.y) * (g0x * texture2D(tex, p0)  +
                        g1x * texture2D(tex, p1)) +
           g1(fuv.y) * (g0x * texture2D(tex, p2)  +
                        g1x * texture2D(tex, p3));
}

