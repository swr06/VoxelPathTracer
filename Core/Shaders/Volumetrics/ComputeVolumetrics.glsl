#version 430 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

// Bayer dither
#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))


layout (location = 0) out vec3 o_Volumetrics;

in vec2 v_TexCoords;

// Participating media 
uniform sampler3D u_VolumetricDensityData;
layout(r8ui, binding = 1) uniform uimage3D u_VolumetricColorData;



uniform sampler2D u_BlueNoise;
uniform sampler2D u_LinearDepthTexture;

uniform vec3 u_ViewerPosition;
uniform vec2 u_Dimensions;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform bool u_Colored;

uniform int u_LightCount;
uniform float u_Time;
uniform float u_Strength;

layout (std430, binding = 2) buffer SSBO_BlockAverageData
{
    vec4 BlockAverageColorData[128];
};

layout (std430, binding = 4) buffer SSBO_LightData
{
    vec4 LightLocations[1024];
};

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

bool IsInVolume(in vec3 pos)
{
    if (pos.x < 0.0f || pos.y < 0.0f || pos.z < 0.0f || 
        pos.x > float(WORLD_SIZE_X - 1) || pos.y > float(WORLD_SIZE_Y - 1) || pos.z > float(WORLD_SIZE_Z - 1))
    {
        return false;    
    }   

    return true;
}

float saturate(float x) {
	return clamp(x, 1e-5f, 1.0f);
}

// Used to create variation
float mod289(float x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 mod289(vec4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 perm(vec4 x){return mod289(((x * 34.0) + 1.0) * x);}

float noise(vec3 p)
{
    vec3 a = floor(p);
    vec3 d = p - a;
    d = d * d * (3.0 - 2.0 * d);
    vec4 b = a.xxyy + vec4(0.0, 1.0, 0.0, 1.0);
    vec4 k1 = perm(b.xyxy);
    vec4 k2 = perm(k1.xyxy + b.zzww);
    vec4 c = k2 + a.zzzz;
    vec4 k3 = perm(c);
    vec4 k4 = perm(c + 1.0);
    vec4 o1 = fract(k3 * (1.0 / 41.0));
    vec4 o2 = fract(k4 * (1.0 / 41.0));
    vec4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    vec2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);
    return o4.y * d.y + o4.x * (1.0 - d.y);
}

// Random rotation matrices 
const mat3 rot1 = mat3(-0.37, 0.36, 0.85,-0.14,-0.93, 0.34,0.92, 0.01,0.4);
const mat3 rot2 = mat3(-0.55,-0.39, 0.74, 0.33,-0.91,-0.24,0.77, 0.12,0.63);
const mat3 rot3 = mat3(-0.71, 0.52,-0.47,-0.08,-0.72,-0.68,-0.7,-0.45,0.56);

// Random rotation to reduce banding 
float simplex3d_fractal(vec3 m) 
{
	// weight * noise(m * rotation matrix) 
    return   0.5333333 * noise(m*rot1)
			+ 0.2666667 * noise(2.0*m*rot2)
			+ 0.1333333 * noise(4.0*m*rot3)
			+ 0.0666667 * noise(8.0*m);
}

float ManhattanDist(vec3 a, vec3 b) {
	return abs(a.x - b.x) + abs(a.y - b.y) + abs(a.z - b.z);
}

// 3D hash3d2 function
float hash3d2(vec3 p)
{
    p  = fract( p*0.3183099+.1 );
	p *= 17.0;
    return fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}

// 3D precedural noise3d2
float noise3d2( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
	
    // interpolate between hashes of adjacent grid points
    return mix(mix(mix( hash3d2(p+vec3(0,0,0)), 
                        hash3d2(p+vec3(1,0,0)),f.x),
                   mix( hash3d2(p+vec3(0,1,0)), 
                        hash3d2(p+vec3(1,1,0)),f.x),f.y),
               mix(mix( hash3d2(p+vec3(0,0,1)), 
                        hash3d2(p+vec3(1,0,1)),f.x),
                   mix( hash3d2(p+vec3(0,1,1)), 
                        hash3d2(p+vec3(1,1,1)),f.x),f.y),f.z);
}

// 3D noise3d2 layered in several octaves
float layeredNoise(in vec3 x) {
    x += vec3(10.0, 5.0, 6.0);
    return 0.6*noise3d2(x*5.0) + 0.4*noise3d2(x*10.0) + 0.2*noise3d2(x*16.0) - 0.2;
}

vec3 sample_triquadratic_gradient_approx(sampler3D channel, vec3 res, vec3 uv) {
    vec3 q = fract(uv * res);
    vec3 cc = 0.5 / res;
    vec3 ww0 = uv - cc;
    vec3 ww1 = uv + cc;
    float nx = texture(channel, vec3(ww1.x, uv.y, uv.z)).r - texture(channel, vec3(ww0.x, uv.y, uv.z)).r;
    float ny = texture(channel, vec3(uv.x, ww1.y, uv.z)).r - texture(channel, vec3(uv.x, ww0.y, uv.z)).r;
    float nz = texture(channel, vec3(uv.x, uv.y, ww1.z)).r - texture(channel, vec3(uv.x, uv.y, ww0.z)).r;
	return vec3(nx, ny, nz);
}

vec3 SampleVolumetricColor(vec3 UV) {
    ivec3 Texel = ivec3(UV*vec3(384.0f,128.0f,384.0f));
    uint BlockID = imageLoad(u_VolumetricColorData, Texel).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 SampleVolumetricColorTexel(ivec3 Texel) {
    uint BlockID = imageLoad(u_VolumetricColorData, Texel).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec4 sample_triquadratic(sampler3D channel, vec3 res, vec3 uv) {
    vec3 q = fract(uv * res);
    vec3 c = (q*(q - 1.0) + 0.5) / res;
    vec3 w0 = uv - c;
    vec3 w1 = uv + c;
    vec4 s = texture(channel, vec3(w0.x, w0.y, w0.z))
    	   + texture(channel, vec3(w1.x, w0.y, w0.z))
    	   + texture(channel, vec3(w1.x, w1.y, w0.z))
    	   + texture(channel, vec3(w0.x, w1.y, w0.z))
    	   + texture(channel, vec3(w0.x, w1.y, w1.z))
    	   + texture(channel, vec3(w1.x, w1.y, w1.z))
    	   + texture(channel, vec3(w1.x, w0.y, w1.z))
		   + texture(channel, vec3(w0.x, w0.y, w1.z));
	return s / 8.0;
}

vec3 SoftwareTriquadraticVolumetrics(vec3 uv) {
    const vec3 res = vec3(384.0f, 128.0f, 384.0f);
    vec3 q = fract(uv * res);
    vec3 cc = 0.5 / res;
    vec3 ww0 = uv - cc;
    vec3 ww1 = uv + cc;
    float nx = SampleVolumetricColor(vec3(ww1.x, uv.y, uv.z)).r - SampleVolumetricColor(vec3(ww0.x, uv.y, uv.z)).r;
    float ny = SampleVolumetricColor(vec3(uv.x, ww1.y, uv.z)).r - SampleVolumetricColor(vec3(uv.x, ww0.y, uv.z)).r;
    float nz = SampleVolumetricColor(vec3(uv.x, uv.y, ww1.z)).r - SampleVolumetricColor(vec3(uv.x, uv.y, ww0.z)).r;
	return vec3(nx, ny, nz);
}

void main() 
{
	const int STEPS = 100;
	const int HALF_STEPS = STEPS/2;

	// Ray properties
	float BaseLinearDepth = texture(u_LinearDepthTexture, v_TexCoords).x;
	BaseLinearDepth = BaseLinearDepth < 0.0f ? 10000.0f : BaseLinearDepth;

	vec3 RayDirection = GetRayDirectionAt(v_TexCoords);
	RayDirection = normalize(RayDirection);

    vec3 WorldPositionSample = u_ViewerPosition + RayDirection * BaseLinearDepth;

    float Distance = BaseLinearDepth;
    Distance = clamp(Distance, 0.0f, 96.0f);

    int Steps = 64;
    float SigmaS = 0.3f * 0.0625f * 0.05f * 5.0f; 

    float StepSize = Distance / float(Steps);
	vec3 BlueNoise = texture(u_BlueNoise, v_TexCoords * (u_Dimensions / vec2(256.0f))).xyz;
    StepSize *= BlueNoise.x;

    vec3 RayPosition = u_ViewerPosition + RayDirection * 0.2f;

    vec3 VolumetricLighting = vec3(0.0f);

    const bool UsePerlinFBM = false;

    for (int CurrentStep = 0 ; CurrentStep < Steps ; CurrentStep++) {
        
        float AirDensity = 1.0f; // Base density at position

        if (UsePerlinFBM) {
            AirDensity = noise(RayPosition);
        }

        if (AirDensity > 0.001f) {
            float SampleSigmaS = SigmaS * AirDensity; 
            float PointDensity = texture(u_VolumetricDensityData, RayPosition.xyz / vec3(384.0f, 128.0f, 384.0f)).x * 24.0f; // Hardware trilinear
            vec3 ScatteringColor = vec3(1.0f);

            if (u_Colored) {
                // Triquadratic interpolation with cubic spline with an approximated gradient curve
                ScatteringColor = clamp(SoftwareTriquadraticVolumetrics(RayPosition/vec3(384.0f,128.0f,384.0f)), 0.0f, 1.0f); 
               // ScatteringColor = clamp(SampleVolumetricColorTexel(ivec3(RayPosition)), 0.0f, 1.0f); 
            }

            vec3 CurrentDensity = vec3(pow(2.0f, 7.0f) * pow(PointDensity, 2.0f))*ScatteringColor;
            vec3 S = SampleSigmaS * CurrentDensity;
            float Transmittance = exp(-(SigmaS*AirDensity) * StepSize); // standard
			vec3 IntegratedScattering = (S - S * Transmittance) / (SigmaS*AirDensity);  // S-St/Se
            VolumetricLighting += IntegratedScattering;
        }

        RayPosition += RayDirection * StepSize;
    }

    const float IsotropicScatter = 0.25 / PI;
    const float ScattersM = pow(2.0f, 4.0f);
    VolumetricLighting *= ScattersM * u_Strength;
    o_Volumetrics = VolumetricLighting;
}