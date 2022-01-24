#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_BloomMips[5];
uniform sampler2D u_BloomBrightTexture;

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

void main() {

	vec3 Bloom[5] = vec3[](vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f));

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

	float Weights[5] = float[5](5.75f, 4.95f, 4.9f, 4.8f, 4.75f);
	const float DetailWeight = 7.2f;

	TotalBloom = (BaseBrightTex * DetailWeight * 1.0f) + TotalBloom;
	TotalBloom = (BaseBrightTex * DetailWeight * 1.0f) + TotalBloom;
	TotalBloom = (pow(Bloom[0], vec3(1.0f / 1.1f)) * Weights[0]) + TotalBloom;
	TotalBloom = (pow(Bloom[1], vec3(1.0f / 1.1f)) * Weights[1]) + TotalBloom;
	TotalBloom = (pow(Bloom[2], vec3(1.0f / 1.05f)) * Weights[2]) + TotalBloom;
	TotalBloom = (pow(Bloom[3], vec3(1.0f / 1.05f)) * Weights[3]) + TotalBloom;
	TotalBloom = (pow(Bloom[4], vec3(1.0f / 1.05f)) * Weights[4]) + TotalBloom;

	float TotalWeights = DetailWeight + Weights[0] + Weights[1] + Weights[2] + Weights[3] + Weights[4];
	TotalBloom /= TotalWeights;


	o_Color = TotalBloom;
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
