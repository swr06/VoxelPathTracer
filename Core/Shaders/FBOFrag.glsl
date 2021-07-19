#version 330 core

layout(location = 0) out vec3 o_Color;
in vec2 v_TexCoords;

uniform sampler2D u_FramebufferTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_BlockIDTex;
uniform bool u_BrutalFXAA;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SamplePositionAt(vec2 txc)
{
	vec3 O = u_InverseView[3].xyz;
	float Dist = texture(u_PositionTexture, txc).r;
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

float quality[12] = float[12] (1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 2.0, 2.0, 2.0, 2.0, 4.0, 8.0);

bool DetectEdge()
{
	vec3 BasePosition = SamplePositionAt(v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec3 BaseColor = texture(u_FramebufferTexture, v_TexCoords).xyz;
	vec2 TexelSize = 1.0f / textureSize(u_FramebufferTexture, 0);
	int BaseBlock = GetBlockAt(v_TexCoords);

	for (int x = -1 ; x <= 1 ; x++)
	{
		for (int y = -2 ; y <= 2 ; y++)
		{
			vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;
			vec3 SamplePosition = SamplePositionAt(SampleCoord).xyz;
			vec3 SampleNormal = texture(u_NormalTexture, SampleCoord).xyz;
			float PositionError = distance(BasePosition, SamplePosition);
			int SampleBlock = GetBlockAt(SampleCoord);

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
	float subpixelQuality = u_BrutalFXAA ? (DetectEdge() ? 3.0 : 1.0) : 0.8f; 
	int iterations = 12;
	vec2 texCoord = v_TexCoords;
	
	vec2 view = 1.0 / vec2(textureSize(u_FramebufferTexture, 0));
	
	float lumaCenter = GetLuminance(color);
	float lumaDown  = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2( 0.0, -1.0) * view, 0.0).rgb);
	float lumaUp    = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2( 0.0,  1.0) * view, 0.0).rgb);
	float lumaLeft  = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2(-1.0,  0.0) * view, 0.0).rgb);
	float lumaRight = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2( 1.0,  0.0) * view, 0.0).rgb);
	
	float lumaMin = min(lumaCenter, min(min(lumaDown, lumaUp), min(lumaLeft, lumaRight)));
	float lumaMax = max(lumaCenter, max(max(lumaDown, lumaUp), max(lumaLeft, lumaRight)));
	
	float lumaRange = lumaMax - lumaMin;
	
	if (lumaRange > max(edgeThresholdMin, lumaMax * edgeThresholdMax)) {
		float lumaDownLeft  = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2(-1.0, -1.0) * view, 0.0).rgb);
		float lumaUpRight   = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2( 1.0,  1.0) * view, 0.0).rgb);
		float lumaUpLeft    = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2(-1.0,  1.0) * view, 0.0).rgb);
		float lumaDownRight = GetLuminance(textureLod(u_FramebufferTexture, texCoord + vec2( 1.0, -1.0) * view, 0.0).rgb);
		
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

		float lumaEnd1 = GetLuminance(textureLod(u_FramebufferTexture, uv1, 0.0).rgb);
		float lumaEnd2 = GetLuminance(textureLod(u_FramebufferTexture, uv2, 0.0).rgb);
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
					lumaEnd1 = GetLuminance(textureLod(u_FramebufferTexture, uv1, 0.0).rgb);
					lumaEnd1 = lumaEnd1 - lumaLocalAverage;
				}
				if (!reached2) {
					lumaEnd2 = GetLuminance(textureLod(u_FramebufferTexture, uv2, 0.0).rgb);
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

	vec3 BasePosition = SamplePositionAt(v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;
	vec3 BaseColor = texture(u_FramebufferTexture, v_TexCoords).xyz;
	float BaseLuma = GetLuminance(BaseColor);
	int BaseBlock = GetBlockAt(v_TexCoords);
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
			vec2 SampleCoord = v_TexCoords + vec2(x,y) * TexelSize;
			vec3 SampleColor = texture(u_FramebufferTexture, SampleCoord).xyz;
			AverageLuminosity += GetLuminance(SampleColor) * float(x != 0 && y != 0) * CurrentWeight;
			BlurredColor += SampleColor * CurrentWeight;
			TotalWeight += 1.0f * CurrentWeight;

			if (x < -1 || x > 1 || y < -1 || y > 1) {
				continue;
			}

			// Edge detection : 

			vec3 SamplePosition = SamplePositionAt(SampleCoord).xyz;
			vec3 SampleNormal = texture(u_NormalTexture, SampleCoord).xyz;
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


void main()
{
	vec3 BaseSample = texture(u_FramebufferTexture, v_TexCoords).rgb;
    vec3 Color = BaseSample;
	const bool CustomFXAA = false;

	if (CustomFXAA) {
		Color = GetFXAACustom();
	} else {
		FXAA311(Color);
	}

	o_Color = pow(Color, vec3(1.0f / 2.2f)); // Gamma correction
}