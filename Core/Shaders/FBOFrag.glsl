#version 330 core
#define INF 100000.0f

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

layout(location = 0) out vec3 o_Color;

uniform sampler2D u_FramebufferTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_BlockIDTex;
uniform sampler2D u_ColorLUT;

uniform bool u_BrutalFXAA;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform vec2 u_Dimensions;

uniform int u_Padding;
uniform int u_SelectedLUT;

uniform bool u_CAS;
uniform bool u_FXAA;

uniform bool u_ColorGrading;
uniform bool u_ColorDither;

uniform bool u_RenderItemCube;

uniform bool u_ExponentiallyMagnifyColorDifferences;

uniform float u_Time;

vec2 TexCoords;

struct HitRecord
{
    float HitDistance;
};

float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}


mat3 RotationMatrix(vec3 axis, float ang)
{
    axis = normalize(axis);
    float s = sin(ang);
    float c = cos(ang);
    float oc = 1.0 - c;
    return mat3(oc*axis.x*axis.x+c,oc*axis.x*axis.y-axis.z*s,oc*axis.z*axis.x+axis.y*s,
                oc*axis.x*axis.y+axis.z*s,oc*axis.y*axis.y+c,oc*axis.y*axis.z-axis.x*s,
                oc*axis.z*axis.x-axis.y*s,oc*axis.y*axis.z+axis.x*s,oc*axis.z*axis.z+c);
}

float RayBoxIntersectionTest(vec3 raypos, vec3 raydir, vec3 boxmin, vec3 boxmax)
{
    float t1 = (boxmin.x - raypos.x) / raydir.x;
    float t2 = (boxmax.x - raypos.x) / raydir.x;
    float t3 = (boxmin.y - raypos.y) / raydir.y;
    float t4 = (boxmax.y - raypos.y) / raydir.y;
    float t5 = (boxmin.z - raypos.z) / raydir.z;
    float t6 = (boxmax.z - raypos.z) / raydir.z;

    float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    if (tmax < 0.0) 
    {
        return INF;
    }

    if (tmin > tmax)
    {
        return INF;
    }

    return tmin;
}

HitRecord IntersectItemCube(vec3 r0, vec3 rD, vec3 BoxMinCoord, vec3 BoxMaxCoord)
{
    HitRecord result;
    result.HitDistance = RayBoxIntersectionTest(r0, rD, BoxMinCoord, BoxMaxCoord);

    float Transversal = step(result.HitDistance, INF);
    
    return result;
}

bool IntersectItemCube() 
{
	if (!u_RenderItemCube) {
		return false;
	}

	vec2 LocalUV;
        LocalUV.x = remap(TexCoords.x, 0.75f, 1.0f, 0.0f, 1.0f);
        LocalUV.y = remap(TexCoords.y, 0.65f, 1.0f, 0.0f, 1.0f);

	if (!(TexCoords.x > 0.8f && TexCoords.y > 0.78f)) {
        return false;
    }

    vec2 NDC = LocalUV * 2.0f - 1.0f;
    vec2 Clip = NDC;
    float SinTime = sin(u_Time);
    mat3 RotationMatrix = RotationMatrix(vec3(1.1f, 3.0f, 1.1f), u_Time * 0.75f);
    vec3 RayDirection = RotationMatrix * vec3(NDC, 1.0f);
    vec3 RayOrigin = RotationMatrix * vec3(0.0f, 1.0f, -1.75f);
    RayDirection += vec3(0.0f, -0.52f, 0.0f);
    RayDirection = normalize(RayDirection);
    HitRecord result = IntersectItemCube(RayOrigin, RayDirection, vec3(0.0f), vec3(1.0f));

    if(result.HitDistance >= INF - 0.025f)
    {
        return false;
    }

    return true;
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

vec3 SampleNormalFromTex(sampler2D samp, vec2 txc) { 
    return GetNormalFromID(texture(samp, txc).x);

}
vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SamplePositionAt(vec2 txc)
{
	vec3 O = u_InverseView[3].xyz;
	float Dist = 1./texture(u_PositionTexture, txc).r;
	return vec4(O + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

int GetBlockAt(vec2 txc)
{
	float id = texture(u_BlockIDTex, txc).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

float GetLuminance(vec3 color) 
{
	return dot(color, vec3(0.299, 0.587, 0.114));
}

float linearToSrgb(float linear){
    float SRGBLo = linear * 12.92;
    float SRGBHi = (pow(abs(linear), 1.0/2.4) * 1.055) - 0.055;
    float SRGB = mix(SRGBHi, SRGBLo, step(linear, 0.0031308));
    return SRGB;
}

float srgbToLinear(float color) {
    float linearRGBLo = color / 12.92;
    float linearRGBHi = pow((color + 0.055) / 1.055, 2.4);
    float linearRGB = mix(linearRGBHi, linearRGBLo, step(color, 0.04045));
    return linearRGB;
}

vec3 linearToSrgb(vec3 linear) {
    vec3 SRGBLo = linear * 12.92;
    vec3 SRGBHi = (pow(abs(linear), vec3(1.0/2.4)) * 1.055) - 0.055;
    vec3 SRGB = mix(SRGBHi, SRGBLo, step(linear, vec3(0.0031308)));
    return SRGB;
}

vec3 srgbToLinear(vec3 color) {
    vec3 linearRGBLo = color / 12.92;
    vec3 linearRGBHi = pow((color + 0.055) / 1.055, vec3(2.4));
    vec3 linearRGB = mix(linearRGBHi, linearRGBLo, step(color, vec3(0.04045)));
    return linearRGB;
}

vec3 Reinhard(vec3 RGB )
{
    return vec3(RGB) / (vec3(1.0f) + GetLuminance(RGB));
}

vec3 InverseReinhard(vec3 RGB)
{
    return RGB / (vec3(1.0f) - GetLuminance(RGB));
}

float GetLuminosityWeightFXAA(vec3 color, bool edge, bool skyedge, vec2 txc) 
{
	if (u_ExponentiallyMagnifyColorDifferences) 
	{
		//color = exp(color*25.0f);
		color = pow(exp(color * 2.4f), vec3(4.0f));
		
	}

	float LuminanceRaw = dot(color, vec3(0.299, 0.587, 0.114));
	return LuminanceRaw;
}

float GetLuminosityWeightFXAANoBias(vec3 color, bool edge, vec2 txc) 
{
	float LuminanceRaw = dot(color, vec3(0.299, 0.587, 0.114));
	return LuminanceRaw;
}

float quality[12] = float[12] (1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 2.0, 2.0, 2.0, 2.0, 4.0, 8.0);

bool DetectEdge(out bool Skyedge)
{
	vec4 BasePosition = SamplePositionAt(TexCoords).xyzw;
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, TexCoords).xyz;
	vec3 BaseColor = texture(u_FramebufferTexture, TexCoords).xyz;
	vec2 TexelSize = 1.0f / textureSize(u_FramebufferTexture, 0);
	int BaseBlock = GetBlockAt(TexCoords);
	Skyedge = false;

	for (int x = -1 ; x <= 1 ; x++)
	{
		for (int y = -1 ; y <= 1 ; y++)
		{
			vec2 SampleCoord = TexCoords + vec2(x, y) * TexelSize;
			vec4 SamplePosition = SamplePositionAt(SampleCoord).xyzw;
			vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
			float PositionError = abs(SamplePosition.w - BasePosition.w);//distance(BasePosition, SamplePosition.xyz);
			int SampleBlock = GetBlockAt(SampleCoord);

			if (SamplePosition.w < 0.0f) {
				Skyedge = true;
			}

			if (BaseNormal != SampleNormal ||
				PositionError > 0.9f ||
				SampleBlock != BaseBlock) 
			{
				return true;
			}
		}
	}

	return false;
}



void FXAA311(inout vec3 color) 
{
	float edgeThresholdMin = 0.03125;
	float edgeThresholdMax = 0.125;
	bool Skyedge = false;
	bool IsAtEdge = DetectEdge(Skyedge);
	float subpixelQuality = IsAtEdge ? 0.935f : 0.0925f; 

	//if (IsAtEdge) {
	//	color = vec3(1.,0.,0.);
	//	return;
	//}

	if (!IsAtEdge) {
		if (IntersectItemCube()) {
			if (u_ExponentiallyMagnifyColorDifferences) {
				return;
			}

			subpixelQuality = 0.825f;
		}
	}

	int iterations = 12;
	vec2 texCoord = TexCoords;

	//if (IsAtEdge) {
	//	color = vec3(1.0f, 0.0f, 0.0f);
	//	return;
	//}
	
	vec2 view = 1.0 / vec2(textureSize(u_FramebufferTexture, 0));
	
	float lumaCenter = GetLuminosityWeightFXAA(color, IsAtEdge, Skyedge, texCoord);
	float lumaDown  = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2( 0.0, -1.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2( 0.0, -1.0) * view);
	float lumaUp    = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2( 0.0,  1.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2( 0.0,  1.0) * view);
	float lumaLeft  = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2(-1.0,  0.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2(-1.0,  0.0) * view);
	float lumaRight = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2( 1.0,  0.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2( 1.0,  0.0) * view);
	
	float lumaMin = min(lumaCenter, min(min(lumaDown, lumaUp), min(lumaLeft, lumaRight)));
	float lumaMax = max(lumaCenter, max(max(lumaDown, lumaUp), max(lumaLeft, lumaRight)));
	
	float lumaRange = lumaMax - lumaMin;
	
	if (lumaRange > max(edgeThresholdMin, lumaMax * edgeThresholdMax)) {
		float lumaDownLeft  = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2(-1.0, -1.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2(-1.0, -1.0) * view);
		float lumaUpRight   = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2( 1.0,  1.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2( 1.0,  1.0) * view);
		float lumaUpLeft    = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2(-1.0,  1.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2(-1.0,  1.0) * view);
		float lumaDownRight = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, texCoord + vec2( 1.0, -1.0) * view, 0.0).rgb, IsAtEdge, Skyedge, texCoord + vec2( 1.0, -1.0) * view);
		
		float lumaDownUp    = lumaDown + lumaUp;
		float lumaLeftRight = lumaLeft + lumaRight;
		
		float lumaLeftCorners  = lumaDownLeft  + lumaUpLeft;
		float lumaDownCorners  = lumaDownLeft  + lumaDownRight;
		float lumaRightCorners = lumaDownRight + lumaUpRight;
		float lumaUpCorners    = lumaUpRight   + lumaUpLeft;
		
		float edgeHorizontal = abs(-2.0 * lumaLeft   + lumaLeftCorners ) +
							   abs(-2.0 * lumaCenter + lumaDownUp      ) * 2.0 +
							   abs(-2.0 * lumaRight  + lumaRightCorners);
		float edgeVertical   = abs(-2.0 * lumaUp     + lumaUpCorners   ) +
							   abs(-2.0 * lumaCenter + lumaLeftRight   ) * 2.0 +
							   abs(-2.0 * lumaDown   + lumaDownCorners );
		
		bool isHorizontal = (edgeHorizontal >= edgeVertical);		
		
		float luma1 = isHorizontal ? lumaDown : lumaLeft;
		float luma2 = isHorizontal ? lumaUp : lumaRight;
		float gradient1 = luma1 - lumaCenter;
		float gradient2 = luma2 - lumaCenter;
		
		bool is1Steepest = abs(gradient1) >= abs(gradient2);
		float gradientScaled = 0.25 * max(abs(gradient1), abs(gradient2));
		
		float stepLength = isHorizontal ? view.y : view.x;

		float lumaLocalAverage = 0.0;

		if (is1Steepest) {
			stepLength = - stepLength;
			lumaLocalAverage = 0.5 * (luma1 + lumaCenter);
		} else {
			lumaLocalAverage = 0.5 * (luma2 + lumaCenter);
		}
		
		vec2 currentUv = texCoord;
		if (isHorizontal) {
			currentUv.y += stepLength * 0.5;
		} else {
			currentUv.x += stepLength * 0.5;
		}
		
		vec2 offset = isHorizontal ? vec2(view.x, 0.0) : vec2(0.0, view.y);
		
		vec2 uv1 = currentUv - offset;
		vec2 uv2 = currentUv + offset;

		float lumaEnd1 = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, uv1, 0.0).rgb, IsAtEdge, Skyedge, uv1);
		float lumaEnd2 = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, uv2, 0.0).rgb, IsAtEdge, Skyedge, uv2);
		lumaEnd1 -= lumaLocalAverage;
		lumaEnd2 -= lumaLocalAverage;
		
		bool reached1 = abs(lumaEnd1) >= gradientScaled;
		bool reached2 = abs(lumaEnd2) >= gradientScaled;
		bool reachedBoth = reached1 && reached2;
		
		if (!reached1) {
			uv1 -= offset;
		}
		if (!reached2) {
			uv2 += offset;
		}
		
		if (!reachedBoth) {
			for(int i = 2; i < iterations; i++) {
				if (!reached1) {
					lumaEnd1 = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, uv1, 0.0).rgb, IsAtEdge, Skyedge, uv1);
					lumaEnd1 = lumaEnd1 - lumaLocalAverage;
				}
				if (!reached2) {
					lumaEnd2 = GetLuminosityWeightFXAA(textureLod(u_FramebufferTexture, uv2, 0.0).rgb, IsAtEdge, Skyedge, uv2);
					lumaEnd2 = lumaEnd2 - lumaLocalAverage;
				}
				
				reached1 = abs(lumaEnd1) >= gradientScaled;
				reached2 = abs(lumaEnd2) >= gradientScaled;
				reachedBoth = reached1 && reached2;

				if (!reached1) {
					uv1 -= offset * quality[i];
				}
				if (!reached2) {
					uv2 += offset * quality[i];
				}
				
				if (reachedBoth) break;
			}
		}
		
		float distance1 = isHorizontal ? (texCoord.x - uv1.x) : (texCoord.y - uv1.y);
		float distance2 = isHorizontal ? (uv2.x - texCoord.x) : (uv2.y - texCoord.y);

		bool isDirection1 = distance1 < distance2;
		float distanceFinal = min(distance1, distance2);

		float edgeThickness = (distance1 + distance2);

		float pixelOffset = - distanceFinal / edgeThickness + 0.5f;
		
		bool isLumaCenterSmaller = lumaCenter < lumaLocalAverage;

		bool correctVariation = ((isDirection1 ? lumaEnd1 : lumaEnd2) < 0.0) != isLumaCenterSmaller;

		float finalOffset = correctVariation ? pixelOffset : 0.0;
		
		float lumaAverage = (1.0 / 12.0) * (2.0 * (lumaDownUp + lumaLeftRight) + lumaLeftCorners + lumaRightCorners);
		float subPixelOffset1 = clamp(abs(lumaAverage - lumaCenter) / lumaRange, 0.0, 1.0);
		float subPixelOffset2 = (-2.0 * subPixelOffset1 + 3.0) * subPixelOffset1 * subPixelOffset1;
		float subPixelOffsetFinal = subPixelOffset2 * subPixelOffset2 * subpixelQuality;

		finalOffset = max(finalOffset, subPixelOffsetFinal);
		
		
		// Compute the final UV coordinates.
		vec2 finalUv = texCoord;
		if (isHorizontal) {
			finalUv.y += finalOffset * stepLength;
		} else {
			finalUv.x += finalOffset * stepLength;
		}

		color = textureLod(u_FramebufferTexture, finalUv, 0.0).rgb;
	}
}

//
// W I P! 
//
vec3 GetFXAACustom()
{
	vec3 BlurredColor = vec3(0.0f);
	vec2 TexelSize = 1.0f / textureSize(u_FramebufferTexture, 0);

	vec3 BasePosition = SamplePositionAt(TexCoords).xyz;
	vec3 BaseNormal = SampleNormalFromTex(u_NormalTexture, TexCoords).xyz;
	vec3 BaseColor = texture(u_FramebufferTexture, TexCoords).xyz;
	float BaseLuma = GetLuminance(BaseColor);
	int BaseBlock = GetBlockAt(TexCoords);
	float AverageLuminosity = 0.0f;
	int AliasedSamples = 0;
	int SamplesTaken = 0;
	bool IsAtEdge = false;

	float TotalWeight = 0.0f;
	float Step = 1.0f;

	for (int xx = -2 ; xx <= 2 ; xx++)
	{
		for (int yy = -2 ; yy <= 2; yy++)
		{	
			float x = xx * Step; 
			float y = yy * Step;
			float CurrentWeight = 1.0f;
			vec2 SampleCoord = TexCoords + vec2(x,y) * TexelSize;
			vec3 SampleColor = texture(u_FramebufferTexture, SampleCoord).xyz;
			AverageLuminosity += GetLuminance(SampleColor) * float(x != 0 && y != 0) * CurrentWeight;
			BlurredColor += SampleColor * CurrentWeight;
			TotalWeight += 1.0f * CurrentWeight;

			if (x < -1 || x > 1 || y < -1 || y > 1) {
				continue;
			}

			// Edge detection : 

			vec3 SamplePosition = SamplePositionAt(SampleCoord).xyz;
			vec3 SampleNormal = SampleNormalFromTex(u_NormalTexture, SampleCoord).xyz;
			float PositionError = distance(BasePosition, SamplePosition);
			float SampleLuma = GetLuminance(SampleColor);
			float LumaDifference = abs(SampleLuma - BaseLuma);
			int SampleBlock = GetBlockAt(SampleCoord);

			LumaDifference = 0.0f;

			if (BaseNormal != SampleNormal ||
				PositionError > 0.9f ||
				LumaDifference >= 0.3f || SampleBlock != BaseBlock) 
			{
				AliasedSamples++;
				IsAtEdge = true;
			}

			SamplesTaken += 1;
		}
	}

	BlurredColor /= TotalWeight;
	AverageLuminosity /= TotalWeight;


	float Mix = float(AliasedSamples) / float(SamplesTaken);
	Mix = clamp(Mix, 0.75f, 1.0f);

	if (IsAtEdge)
	{
		return mix(BaseColor, BlurredColor, Mix);
	}

	return BaseColor;

}

vec3 Lookup(vec3 color) 
{
    const vec2 InverseSize = vec2(1.0f / 512.0f, 1.0f / 5120.0f);
    mat2 Grid = mat2(vec2(1.0f, InverseSize.y * 512), vec2(0.0f, u_SelectedLUT * InverseSize.y * 512)); 

    // calculate blue component -> used to figure out the quad
    float BlueComponent = color.b * 63.0f;

    // Get Quad 
    vec4 Quad = vec4(0.0f);
    Quad.y = floor(floor(BlueComponent) * 0.125f);
    Quad.x = floor(BlueComponent) - (Quad.y * 8.0f);
    Quad.w = floor(ceil(BlueComponent) * 0.125f);
    Quad.z = ceil(BlueComponent) - (Quad.w * 8.0f);

    // calculate sample pos
    vec4 SamplePosition = (Quad * 0.125f) + (0.123046875f * color.rg).xyxy + 0.0009765625f;

    // fetch
    vec3 Fetch1 = texture2D(u_ColorLUT, SamplePosition.xy * Grid[0] + Grid[1]).rgb;
    vec3 Fetch2 = texture2D(u_ColorLUT, SamplePosition.zw * Grid[0] + Grid[1]).rgb;
    
    // mix
    return mix(Fetch1, Fetch2, fract(BlueComponent));
}

void BasicColorDither(inout vec3 color)
{
    color += bayer128(gl_FragCoord.xy) / 128.0f;
	const vec2 LestynCoefficients = vec2(171.0f, 231.0f);
    vec3 Lestyn = vec3(dot(LestynCoefficients, gl_FragCoord.xy));
    Lestyn = fract(Lestyn.rgb / vec3(103.0f, 71.0f, 97.0f));
    color += Lestyn.rgb / 255.0f;
}

void main()
{

	// clip the screen 
	TexCoords = vec2(gl_FragCoord.xy);
	TexCoords += u_Padding/2;
	TexCoords = TexCoords / vec2(u_Dimensions + float(u_Padding));

	vec3 BaseSample = texture(u_FramebufferTexture, TexCoords).rgb;
	vec3 ViewerPos = u_InverseView[3].xyz;
	vec3 BasePos = SamplePositionAt(TexCoords).xyz;

    vec3 Color = BaseSample;
	const bool CustomFXAA = false;
	bool fxaa = false;

	if (u_FXAA) {
		FXAA311(Color);
	}

	o_Color = Color;

	if (!u_CAS) {
		o_Color = linearToSrgb(Color); // Gamma correction
		o_Color = clamp(o_Color, 0.0f, 1.0f);

		if (u_ColorGrading) {
			o_Color = Lookup(o_Color);
		}

		if (u_ColorDither) {
			BasicColorDither(o_Color);
		}
	}

	else {
		o_Color = clamp(o_Color, 0.0f, 1.0f);
	}
}