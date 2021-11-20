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

vec4 bicubicTexture(sampler2D tex, vec2 coord, float res) {
  vec2 resolutionScale = resolution * res;

  coord *= resolutionScale;

  float fx = fract(coord.x);
  float fy = fract(coord.y);
  coord.x -= fx;
  coord.y -= fy;

  vec4 xcubic = cubic(fx);
  vec4 ycubic = cubic(fy);

  vec4 c = vec4(coord.x - 0.5, coord.x + 1.5, coord.y - 0.5, coord.y + 1.5);
  vec4 s = vec4(xcubic.x + xcubic.y, xcubic.z + xcubic.w, ycubic.x + ycubic.y, ycubic.z + ycubic.w);
  vec4 offset = c + vec4(xcubic.y, xcubic.w, ycubic.y, ycubic.w) / s;

  vec4 sample0 = texture2DLod(tex, vec2(offset.x, offset.z) / resolutionScale, 0.0);
  vec4 sample1 = texture2DLod(tex, vec2(offset.y, offset.z) / resolutionScale, 0.0);
  vec4 sample2 = texture2DLod(tex, vec2(offset.x, offset.w) / resolutionScale, 0.0);
  vec4 sample3 = texture2DLod(tex, vec2(offset.y, offset.w) / resolutionScale, 0.0);

  float sx = s.x / (s.x + s.y);
  float sy = s.z / (s.z + s.w);

  return mix( mix(sample3, sample2, sx), mix(sample1, sample0, sx), sy);
}

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

vec4 textureSmoothGrad(sampler2D t, vec2 x, vec2 textureSize, mat2 derivatives) {
    x *= vec2(textureSize);

    vec2 p = floor(x);
    vec2 f = fract(x);

    vec4 a = texture2DGrad(t, (p                 ) / vec2(textureSize), derivatives[0], derivatives[1]);
    vec4 b = texture2DGrad(t, (p + vec2(1.0, 0.0)) / vec2(textureSize), derivatives[0], derivatives[1]);
    vec4 c = texture2DGrad(t, (p + vec2(0.0, 1.0)) / vec2(textureSize), derivatives[0], derivatives[1]);
    vec4 d = texture2DGrad(t, (p + vec2(1.0, 1.0)) / vec2(textureSize), derivatives[0], derivatives[1]);

    return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}