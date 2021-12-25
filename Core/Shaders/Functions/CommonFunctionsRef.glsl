#version 430 core

vec4 softwareBilinear(sampler2D tex, vec2 uv)
{
    vec2 texSize = textureSize(tex, 0).xy;
	vec2 pos = uv * texSize - 0.5;
    vec2 f = fract(pos);
    vec2 pos_top_left = floor(pos);
    vec4 tl = texture(tex, (pos_top_left + vec2(0.5, 0.5)) / texSize, -100.0);
    vec4 tr = texture(tex, (pos_top_left + vec2(1.5, 0.5)) / texSize, -100.0);
    vec4 bl = texture(tex, (pos_top_left + vec2(0.5, 1.5)) / texSize, -100.0);
    vec4 br = texture(tex, (pos_top_left + vec2(1.5, 1.5)) / texSize, -100.0);
    vec4 ret = mix(mix(tl, tr, f.x), mix(bl, br, f.x), f.y);
    return ret;
}

vec4 nearest(sampler2D tex, vec2 uv) {
    vec2 res = textureSize(tex,0).xy;
    uv *= res;
    uv = floor(uv)+0.5;
    uv /= res;
    return textureLod(tex, uv, 0.0);
}

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

float log10(float x){
    return log(x) / log(10.0);
}
vec2 log10(vec2 x){
    return log(x) / log(10.0);
}
vec3 log10(vec3 x){
    return log(x) / log(10.0);
}

float max3(float x, float y, float z){
    return max(max(x,y),z);
}
float min3(float x, float y, float z){
    return min(min(x,y),z);
}
float max3(vec3 a){
    return max(max(a.x,a.y),a.z);
}
float min3(vec3 a){
    return min(min(a.x,a.y),a.z);
}
float max0(float x){
    return max(x, 0);
}
float clamp01(float x){
    return clamp(x, 0., 1.);
}
vec2 clamp01(vec2 x){
    return clamp(x, 0., 1.);
}
vec3 clamp01(vec3 x){
    return clamp(x, 0., 1.);
}

vec3 toSRGB(vec3 color) {
	return mix(color * 12.92, 1.055 * pow(color, vec3(1.0 / 2.4)) - 0.055, vec3(greaterThan(color, vec3(0.0031308))));
}

vec3 toLinear(vec3 color) {
	return mix(color / 12.92, pow((color + 0.055) / 1.055, vec3(2.4)), vec3(greaterThan(color, vec3(0.04045))));
}

float Luminance(in vec3 color)
{
	return dot(color.rgb, vec3(0.2125f, 0.7154f, 0.0721f));
}

float LuminanceAccurate(in vec3 color) {
    return dot(color, vec3(0.2722287168, 0.6740817658, 0.0536895174));
}

#define rcp(x) (1.0/x)
#define expf(x) exp2((x) * RCPLOG2)
#define pow2(x) ((x) * (x))
#define pow3(x) (pow2(x) * x)
#define pow4(x) (pow2(x) * pow2(x))
#define pow5(x) (pow2(x) * pow2(x)*x)
#define pow6(x) (pow2(x) * pow2(x) * pow2(x))
#define pow8(x) (pow4(x) * pow4(x))
#define pow12(x) (pow6(x) * pow6(x))
#define pow16(x) (pow8(x) * pow8(x))
#define pow32(x) (pow16(x) * pow16(x))
#define pow64(x) (pow32(x) * pow32(x))
#define pow128(x) (pow64(x) * pow64(x))
#define pow256(x) (pow128(x) * pow128(x))



const float seedDelta = 0.001;

float hash1(inout float seed) 
{
	return fract(sin(seed += seedDelta)*43758.5453123);
}

vec2 hash2(inout float seed) 
{
	return fract(sin(vec2(seed+=seedDelta, seed+=seedDelta))*vec2(43758.5453123, 22578.1459123));
}

vec3 hash3(inout float seed) 
{
	return fract(sin(vec3(seed+=seedDelta, seed+=seedDelta, seed+=seedDelta))*vec3(43758.5453123, 22578.1459123, 19642.3490423));
}

vec4 hash4(inout float seed) 
{
	return fract(sin(vec4(seed+=seedDelta, seed+=seedDelta, seed+=seedDelta, seed+=seedDelta))*vec4(43758.5453123, 22578.1459123, 19642.3490423, 71248.15211892));
}

uint hashi(uint x, uint y, uint z) {
	x += x >> 11;
	x ^= x << 7;
	x += y;
	x ^= x << 3;
	x += z ^ (x >> 14);
	x ^= x << 6;
	x += x >> 15;
	x ^= x << 5;
	x += x >> 12;
	x ^= x << 9;
	return x;
}

float random(vec2 pos, float time) {
	uint mantissaMask = 0x007FFFFFu;
	uint one = 0x3F800000u;
	uvec3 u = floatBitsToUint(vec3(pos, time));
	uint h = hashi(u.x, u.y, u.z);
	return uintBitsToFloat((h & mantissaMask) | one) - 1.0;
}

#define icubeSmooth(x) (x * x) * (3.0 - 2.0 * x)

float cubeSmooth(float x) {
    return icubeSmooth(x);
}
vec2 cubeSmooth(vec2 x) {
    return icubeSmooth(x);
}

float cube_smooth(float x) {
    return (x * x) * (3.0 - 2.0 * x);
}

// by inigo quilez
vec4 smoothfilter(in sampler2D tex, in vec2 uv) {
    vec2 resolution = textureSize(tex, 0);
	uv = uv*resolution + 0.5;
	vec2 iuv = floor(uv);
	vec2 fuv = fract(uv);
	uv = iuv + fuv*fuv*fuv*(fuv*(fuv*6.0-15.0)+10.0);
	uv = (uv - 0.5)/resolution;
	return textureLod(tex, uv, 0.0);
}

// by inigo quilez
vec4 textureSmooth(sampler2D t, vec2 x, vec2 textureSize) {
    x *= vec2(textureSize);

    vec2 p = floor(x);
    vec2 f = fract(x);

    vec4 a = texture2D(t, (p                 ) / vec2(textureSize));
    vec4 b = texture2D(t, (p + vec2(1.0, 0.0)) / vec2(textureSize));
    vec4 c = texture2D(t, (p + vec2(0.0, 1.0)) / vec2(textureSize));
    vec4 d = texture2D(t, (p + vec2(1.0, 1.0)) / vec2(textureSize));

    return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}

// bicubic bspline
vec4 cubic(float x) {
  float x2 = x * x;
  float x3 = x2 * x;
  vec4 w;
  w.x =   -x3 + 3*x2 - 3*x + 1;
  w.y =  3*x3 - 6*x2       + 4;
  w.z = -3*x3 + 3*x2 + 3*x + 1;
  w.w =  x3;
  return w / 6.f;
}

vec4 bicubicTexture(sampler2D tex, vec2 coord, vec2 resolution) {
  coord *= resolution;

  float fx = fract(coord.x);
  float fy = fract(coord.y);
  coord.x -= fx;
  coord.y -= fy;

  vec4 xcubic = cubic(fx);
  vec4 ycubic = cubic(fy);

  vec4 c = vec4(coord.x - 0.5, coord.x + 1.5, coord.y - 0.5, coord.y + 1.5);
  vec4 s = vec4(xcubic.x + xcubic.y, xcubic.z + xcubic.w, ycubic.x + ycubic.y, ycubic.z + ycubic.w);
  vec4 offset = c + vec4(xcubic.y, xcubic.w, ycubic.y, ycubic.w) / s;

  vec4 sample0 = texture2D(tex, vec2(offset.x, offset.z) / resolution);
  vec4 sample1 = texture2D(tex, vec2(offset.y, offset.z) / resolution);
  vec4 sample2 = texture2D(tex, vec2(offset.x, offset.w) / resolution);
  vec4 sample3 = texture2D(tex, vec2(offset.y, offset.w) / resolution);

  float sx = s.x / (s.x + s.y);
  float sy = s.z / (s.z + s.w);

  return mix( mix(sample3, sample2, sx), mix(sample1, sample0, sx), sy);
}

void main() {}