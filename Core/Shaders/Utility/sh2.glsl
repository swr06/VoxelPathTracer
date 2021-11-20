// From q2rtx

struct SphericalHarmonics {
    vec4 Y;
    vec2 CoCg;
};

SphericalHarmonics InitSH() {
    return SphericalHarmonics(vec4(0), vec2(0));
}

void accumulateSH(inout SphericalHarmonics data, SphericalHarmonics b, float scale) {
    data.Y     += b.Y * scale;
    data.CoCg  += b.CoCg * scale;
}
void scaleSH(inout SphericalHarmonics data, float mult) {
    data.Y     *= mult;
    data.CoCg  *= mult;
}
void divideSH(inout SphericalHarmonics data, float denom) {
    data.Y     /= denom;
    data.CoCg  /= denom;
}
SphericalHarmonics mixSH(SphericalHarmonics a, SphericalHarmonics b, float s) {
    return SphericalHarmonics(mix(a.Y, b.Y, vec4(s)), mix(a.CoCg, b.CoCg, vec2(s)));
}


SphericalHarmonics IrradianceToSH(vec3 Color, vec3 Dir) {
    //return SphericalHarmonics(vec4(Color, getLuma(Color)), vec2(0));

    float   Co      = Color.r - Color.b;
    float   t       = Color.b + Co * 0.5;
    float   Cg      = Color.g - t;
    float   Y       = max(t + Cg * 0.5, 0.0);

    float   L00     = 0.282095;
    float   L1_1    = 0.488603 * Dir.y;
    float   L10     = 0.488603 * Dir.z;
    float   L11     = 0.488603 * Dir.x;

    return SphericalHarmonics(vec4(L11, L1_1, L10, L00)*Y, vec2(Co, Cg));
}

vec3 ProjectIrradiance(SphericalHarmonics SH, vec3 N) {
    float D     = dot(SH.Y.xyz, N);
    float Y     = max0(2.0 * (1.023326 * D + 0.886226 * SH.Y.w));

    SH.CoCg    *= Y * 0.282095 / (SH.Y.w + 1e-6);

    float T     = Y - SH.CoCg.y * 0.5;
    float G     = SH.CoCg.y + T;
    float B     = T - SH.CoCg.x * 0.5;
    float R     = B + SH.CoCg.x;

    return max0(vec3(R, G, B));
}
vec3 ToIrradiance(SphericalHarmonics SH) {
    float Y         = SH.Y.w / 0.282095;
    float T         = Y - SH.CoCg.y * 0.5;
    float G         = SH.CoCg.y + T;
    float B         = T - SH.CoCg.x * 0.5;
    float R         = B + SH.CoCg.x;

    return max0(vec3(R, G, B));
}
float IrradianceLuma(SphericalHarmonics SH) {
    return SH.Y.w / 0.282095;
}