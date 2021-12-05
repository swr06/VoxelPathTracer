#version 430 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

// Bayer dithering functions
#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))


layout (location = 0) out vec4 o_Volumetrics;

in vec2 v_TexCoords;

//layout(r8ui, binding = 1) uniform uimage3D u_VolumetricColorData;



uniform sampler2D u_BlueNoise;
uniform sampler2D u_LinearDepthTexture;

// Participating media 
uniform usampler3D u_VolumetricColorDataSampler;
uniform sampler3D u_VolumetricDensityData;

uniform bool u_UsePerlinNoiseForOD;
uniform bool u_PointVolTriquadraticDensityInterpolation;

uniform vec3 u_ViewerPosition;
uniform vec2 u_Dimensions;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform bool u_Colored;
uniform bool u_UseBayer;

uniform int u_LightCount;
uniform float u_Time;
uniform float u_Strength;
uniform bool u_FractalSimplexOD;

uniform bool u_GroundTruthColorInterpolation = false;

// Ssbos 
layout (std430, binding = 2) buffer SSBO_BlockAverageData
{
    vec4 BlockAverageColorData[128]; // Returns the average color per block type 
};

// 
//layout (std430, binding = 4) buffer SSBO_LightData
//{
//    vec4 LightLocations[1024]; 
//};

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


// Noise FBMs 
// Perlin noise
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


// Simplex noise 
vec3 random3(vec3 c) 
{
	float j = 4096.0*sin(dot(c,vec3(17.0, 59.4, 15.0)));
	vec3 r;
	r.z = fract(512.0*j);
	j *= .125;
	r.x = fract(512.0*j);
	j *= .125;
	r.y = fract(512.0*j);
	return r-0.5;
}

// returns simplex noise for a single point 
float simplex3d(vec3 p) {
    const float F3 =  0.3333333f;
    const float G3 =  0.1666667f;
    vec3 s = floor(p + dot(p, vec3(F3)));
    vec3 x = p - s + dot(s, vec3(G3));
    vec3 e = step(vec3(0.0), x - x.yzx);
    vec3 i1 = e*(1.0 - e.zxy);
    vec3 i2 = 1.0 - e.zxy*(1.0 - e);
    vec3 x1 = x - i1 + G3;
    vec3 x2 = x - i2 + 2.0*G3;
    vec3 x3 = x - 1.0 + 3.0*G3;
    vec4 w, d;
    w.x = dot(x, x);
    w.y = dot(x1, x1);
    w.z = dot(x2, x2);
    w.w = dot(x3, x3);
    w = max(0.6 - w, 0.0);
    d.x = dot(random3(s), x);
    d.y = dot(random3(s + i1), x1);
    d.z = dot(random3(s + i2), x2);
    d.w = dot(random3(s + 1.0), x3);
    w *= w;
    w *= w;
    d *= w;
    return dot(d, vec4(52.0));
}

const mat3 RotationMatrix_1 = mat3(-0.37f, 0.36f, 0.85f,-0.14f,-0.93f, 0.34f,0.92f, 0.01f,0.40f);
const mat3 RotationMatrix_2 = mat3(-0.55f,-0.39f, 0.74f, 0.33f,-0.91f,-0.24f,0.77f, 0.12f,0.63f);
const mat3 RotationMatrix_3 = mat3(-0.71f, 0.52f,-0.47f,-0.08f,-0.72f,-0.68f,-0.7f,-0.45f,0.56f);

// Rotate each sample octave to reduce artifacts 
// Works super well usually 
float SimplexFractalNoise(vec3 m) {
    return  4.0f *(  0.5333333f * simplex3d(m * RotationMatrix_1)
			+ 0.2666667f * simplex3d(2.0f * m * RotationMatrix_2)
			+ 0.1333333f * simplex3d(4.0f * m * RotationMatrix_3)
			+ 0.0666667f * simplex3d(8.0f * m));
}


float ScatterIntegral(float OD, float Coefficent)
{
    float InverseLog2 = 1.0f / log(2.0);
    float T = -Coefficent * InverseLog2; 
    float B = -1.0f / Coefficent;
    float C =  1.0f / Coefficent;
    //return exp(T * OD) * B + C;
    return exp2(T * OD) * B + C;
}


// Returns the optical depth for a point using a perlin/simplex fractal FBM
float OpticalDepth(vec3 p)
{
    if (true) {

        //float OD = pow(dot(noise(p - Time),vec3(0.625,0.25,0.625),2.2)*4.0 ;
        //OD -= pow(dot(noise(p - Time*2.0f),vec3(0.625,0.25,0.625),8.2)*4.0 ;
        //OD += pow(dot(noise(p + Time*4.0f),vec3(0.625,0.25,0.625),4.2)*4.0 ;
        //float Integral = ScatterIntegral(OD,1.1f);

        float Time = mod(u_Time * (0.9f / 1.0f), 100.0f);
        float OD = noise(p - Time) * 0.5f; 
        p *= 1.2f;
        OD += noise(p * 2.0f + Time) * 0.25f;
        OD += noise(p * 4.0f - Time) * 0.125f;
        float OpticalDepthSquared = OD * OD;
	    return (OpticalDepthSquared * 4.0f + 0.25f) * (0.9f / 1.0f);
    }

    else {
        float Time = u_Time * 0.9f;
        Time = mod(Time, 100.0f);
        p*=3.0f;
        float OD = SimplexFractalNoise(vec3(p.x - mod(Time, 100.0f), p.y + Time * 0.2f, p.z + Time)) * 1.0f;
	    return (OD*OD * 4.0f + 0.25f) * (0.9f / 1.0f) * (0.9f / 1.0f);
    }
}

float MarchODShadow(vec3 Point, vec3 Lo)
{
    if (!u_UsePerlinNoiseForOD) {
        return 1.0f;
    }

    int Steps = 4;
    float StepSize = 1.0 / float(Steps);
    vec3 Increment = Lo * StepSize;
    vec3 RayPosition = Point;
    float Transmittance = float(1.0);
    
    for (int CurrentStep = 0; CurrentStep < Steps; CurrentStep++)
    {
        float OD_At = OpticalDepth(RayPosition);
		RayPosition += Increment;
        Transmittance += OD_At;
    }
    
    return exp2(-Transmittance * StepSize); //exp(-Transmittance * StepSize * (0.9f / 1.0f))
}


vec3 SampleVolumetricColor(vec3 UV) {
    uint BlockID = texture(u_VolumetricColorDataSampler, UV).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 SampleVolumetricColor(vec3 UV, float D) {
    uint BlockID = texture(u_VolumetricColorDataSampler, UV+D*0.5f*(1.0f/vec3(384.0f,128.0f,384.0f))).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);
}   

vec3 SampleVolumetricColorTexel(ivec3 Texel) {
    uint BlockID = texelFetch(u_VolumetricColorDataSampler, Texel, 0).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 SampleVolumetricColorTexel(ivec3 Texel, int LOD) {
    uint BlockID = texelFetch(u_VolumetricColorDataSampler, Texel, LOD).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 InterpolateVolumetricColorGT(vec3 uv) // Ground truth interpolation
{
    vec3 res = vec3(384.0f, 128.0f, 384.0f);
    vec3 q = fract(uv * res);
    ivec3 t = ivec3(uv * res);
    ivec3 e = ivec3(-1, 0, 1);
    vec3 q0 = (q+1.0)/2.0;
    vec3 q1 = q/2.0;	
    vec3 s000 = SampleVolumetricColorTexel(t + e.xxx, 0);
    vec3 s001 = SampleVolumetricColorTexel(t + e.xxy, 0);
    vec3 s002 = SampleVolumetricColorTexel(t + e.xxz, 0);
    vec3 s012 = SampleVolumetricColorTexel(t + e.xyz, 0);
    vec3 s011 = SampleVolumetricColorTexel(t + e.xyy, 0);
    vec3 s010 = SampleVolumetricColorTexel(t + e.xyx, 0);
    vec3 s020 = SampleVolumetricColorTexel(t + e.xzx, 0);
    vec3 s021 = SampleVolumetricColorTexel(t + e.xzy, 0);
    vec3 s022 = SampleVolumetricColorTexel(t + e.xzz, 0);
    vec3 y00 = mix(mix(s000, s001, q0.z), mix(s001, s002, q1.z), q.z);
    vec3 y01 = mix(mix(s010, s011, q0.z), mix(s011, s012, q1.z), q.z);
    vec3 y02 = mix(mix(s020, s021, q0.z), mix(s021, s022, q1.z), q.z);
	vec3 x0 = mix(mix(y00, y01, q0.y), mix(y01, y02, q1.y), q.y);
    vec3 s122 = SampleVolumetricColorTexel(t + e.yzz, 0);
    vec3 s121 = SampleVolumetricColorTexel(t + e.yzy, 0);
    vec3 s120 = SampleVolumetricColorTexel(t + e.yzx, 0);
    vec3 s110 = SampleVolumetricColorTexel(t + e.yyx, 0);
    vec3 s111 = SampleVolumetricColorTexel(t + e.yyy, 0);
    vec3 s112 = SampleVolumetricColorTexel(t + e.yyz, 0);
    vec3 s102 = SampleVolumetricColorTexel(t + e.yxz, 0);
    vec3 s101 = SampleVolumetricColorTexel(t + e.yxy, 0);
    vec3 s100 = SampleVolumetricColorTexel(t + e.yxx, 0);
    vec3 y10 = mix(mix(s100, s101, q0.z), mix(s101, s102, q1.z), q.z);
    vec3 y11 = mix(mix(s110, s111, q0.z), mix(s111, s112, q1.z), q.z);
    vec3 y12 = mix(mix(s120, s121, q0.z), mix(s121, s122, q1.z), q.z);
    vec3 x1 = mix(mix(y10, y11, q0.y), mix(y11, y12, q1.y), q.y);
    vec3 s200 = SampleVolumetricColorTexel(t + e.zxx, 0);
    vec3 s201 = SampleVolumetricColorTexel(t + e.zxy, 0);
    vec3 s202 = SampleVolumetricColorTexel(t + e.zxz, 0);
    vec3 s212 = SampleVolumetricColorTexel(t + e.zyz, 0);
    vec3 s211 = SampleVolumetricColorTexel(t + e.zyy, 0);
    vec3 s210 = SampleVolumetricColorTexel(t + e.zyx, 0);
    vec3 s220 = SampleVolumetricColorTexel(t + e.zzx, 0);
    vec3 s221 = SampleVolumetricColorTexel(t + e.zzy, 0);
    vec3 s222 = SampleVolumetricColorTexel(t + e.zzz, 0);
    vec3 y20 = mix(mix(s200, s201, q0.z), mix(s201, s202, q1.z), q.z);
    vec3 y21 = mix(mix(s210, s211, q0.z), mix(s211, s212, q1.z), q.z);
    vec3 y22 = mix(mix(s220, s221, q0.z), mix(s221, s222, q1.z), q.z);
    vec3 x2 = mix(mix(y20, y21, q0.y), mix(y21, y22, q1.y), q.y);
    return mix(mix(x0, x1, q0.x), mix(x1, x2, q1.x), q.x);
}

// Custom Triquadratic 3D interpolation 

bool NoInterp = false;
vec3 CustomTriquadraticVolumeInterp(vec3 UV, vec3 BayerNormalized) 
{ 
    if (NoInterp) {
        return SampleVolumetricColor(UV);
    }

    if (u_GroundTruthColorInterpolation) {
        return InterpolateVolumetricColorGT(UV);
    }

    const vec3 Resolution = vec3(384.0f, 128.0f, 384.0f);
    vec3 FractTexel = fract(UV * Resolution);
    vec3 LinearOffset = (FractTexel * (FractTexel - 1.0f) + 0.5f) / Resolution;
    vec3 W0 = UV - LinearOffset;
    vec3 W1 = UV + LinearOffset;

    const float DitherWeights[4] = float[4](1.0f, 0.5f, 0.25f, 0.2f);
    const float GlobalDitherNoiseWeight = 1.414f;

    vec3 Interpolated = SampleVolumetricColor(vec3(W0.x, W0.y, W0.z) + BayerNormalized * DitherWeights[0] * GlobalDitherNoiseWeight)
    	   + SampleVolumetricColor(vec3(W1.x, W0.y, W0.z) - BayerNormalized * DitherWeights[0] * GlobalDitherNoiseWeight)
    	   + SampleVolumetricColor(vec3(W1.x, W1.y, W0.z) + BayerNormalized * DitherWeights[1] * GlobalDitherNoiseWeight)
    	   + SampleVolumetricColor(vec3(W0.x, W1.y, W0.z) - BayerNormalized * DitherWeights[1] * GlobalDitherNoiseWeight)
    	   + SampleVolumetricColor(vec3(W0.x, W1.y, W1.z) + BayerNormalized * DitherWeights[2] * GlobalDitherNoiseWeight)
    	   + SampleVolumetricColor(vec3(W1.x, W1.y, W1.z) - BayerNormalized * DitherWeights[2] * GlobalDitherNoiseWeight)
    	   + SampleVolumetricColor(vec3(W1.x, W0.y, W1.z) + BayerNormalized * DitherWeights[3] * GlobalDitherNoiseWeight)
		   + SampleVolumetricColor(vec3(W0.x, W0.y, W1.z) - BayerNormalized * DitherWeights[3] * GlobalDitherNoiseWeight);

	return max(Interpolated / 8.0, 0.00000001f);
}

float SampleTriQuadraticDensity(vec3 UV) 
{
    vec3 LPVResolution = vec3(384.0f, 128.0f, 384.0f);
    vec3 FractTexel = fract(UV * LPVResolution);
    vec3 LinearOffset = (FractTexel * (FractTexel - 1.0f) + 0.5f) / LPVResolution;
    vec3 W0 = UV - LinearOffset;
    vec3 W1 = UV + LinearOffset;
    float Density = texture(u_VolumetricDensityData, vec3(W0.x, W0.y, W0.z)).x
    	          + texture(u_VolumetricDensityData, vec3(W1.x, W0.y, W0.z)).x
    	          + texture(u_VolumetricDensityData, vec3(W1.x, W1.y, W0.z)).x
    	          + texture(u_VolumetricDensityData, vec3(W0.x, W1.y, W0.z)).x
    	          + texture(u_VolumetricDensityData, vec3(W0.x, W1.y, W1.z)).x
    	          + texture(u_VolumetricDensityData, vec3(W1.x, W1.y, W1.z)).x
    	          + texture(u_VolumetricDensityData, vec3(W1.x, W0.y, W1.z)).x
		          + texture(u_VolumetricDensityData, vec3(W0.x, W0.y, W1.z)).x;
	return max(Density / 8.0, 0.00000001f);
}



void main() 
{
    const int Steps = 64;

	// Ray properties
	float BaseLinearDepth = texture(u_LinearDepthTexture, v_TexCoords).x;
	BaseLinearDepth = BaseLinearDepth < 0.0f ? 10000.0f : BaseLinearDepth;

	vec3 RayDirection = GetRayDirectionAt(v_TexCoords);
	RayDirection = normalize(RayDirection);

    vec3 WorldPositionSample = u_ViewerPosition + RayDirection * BaseLinearDepth;

    float Distance = BaseLinearDepth;
    Distance = clamp(Distance, 0.0f, 96.0f);

    // 
    float SigmaS = 0.00468f; 
    const float SigmaA = 0.1f;
    const float SigmaT = 1.25f;

    float StepSize = Distance / float(Steps);
	vec3 DitherNoise = vec3(1.0f);
    
    // Bayer is much more noiseless
    // Easier to denoise as well!

    if (u_UseBayer) {
        DitherNoise = vec3(bayer256(gl_FragCoord.xy));
        DitherNoise.y = bayer128(gl_FragCoord.xy);
    }

    else {
        DitherNoise = texture(u_BlueNoise, v_TexCoords * (u_Dimensions / vec2(256.0f))).xyz;
    }

    StepSize *= DitherNoise.x;

    vec3 RayPositionActual = u_ViewerPosition + RayDirection * 0.2f;

    vec3 VolumetricLighting = vec3(0.0f);

    bool UsePerlinFBM = u_UsePerlinNoiseForOD;
    float OutputTransmittance = 1.0f;
    float OutputDensitySum = 1.0f;

    vec3 SampleDither = vec3(bayer128(gl_FragCoord.xy), bayer16(gl_FragCoord.xy)*0.4, bayer16(gl_FragCoord.xy)*0.7);
    vec3 VolumetricColorDither = vec3(bayer128(gl_FragCoord.xy));
    VolumetricColorDither *= 1.0f/vec3(384.0f,128.0f,384.0f);

    //float CurrentTransmittance = 1.0f;
    //float TotalTransmittance = 1.0f;

    float pdM = mix(24.0f, 24.0f+4.5f, float(u_PointVolTriquadraticDensityInterpolation));

    for (int CurrentStep = 0 ; CurrentStep < Steps ; CurrentStep++) {
        
        float AirDensity = 1.0f; // Base density at position

        vec3 RayPosition = RayPositionActual + vec3(0.0f, -0.225f, 0.0f);

        if (UsePerlinFBM) {
            AirDensity = clamp(pow(OpticalDepth(RayPosition*0.95f),1.52525f)/1.4144f,0.0f,1.0f);
        }

        if (AirDensity > 0.001f) {

            float SampleSigmaS = SigmaS * AirDensity; 
            vec3 DensitySamplePosition = (vec3(RayPosition.xyz-vec3(0.125f,0.0f,0.125f))*(1.0f/vec3(384.0f,128.0f,384.0f)))+(SampleDither*0.125f*(1.0f/vec3(384.0f,128.0f,384.0f)));
            float PointDensity;

            if (u_PointVolTriquadraticDensityInterpolation) {
                PointDensity = SampleTriQuadraticDensity(DensitySamplePosition).x; // Triquadratic + Trilinear
            }

            else {
                PointDensity = texture(u_VolumetricDensityData, DensitySamplePosition).x; // Hardware trilinear
            }

            PointDensity *= pdM;

            OutputDensitySum += PointDensity;
            vec3 ScatteringColor = vec3(1.0f);

            if (u_Colored) {
                // Triquadratic interpolation 
                vec3 Normalized = RayPosition * (1.0f/vec3(384.0f,128.0f,384.0f));
               // Normalized += SampleDither * 0.450f * (1.0f/vec3(384.0f,128.0f,384.0f)); // Jitter ray sample position
                ScatteringColor = clamp(CustomTriquadraticVolumeInterp(Normalized,VolumetricColorDither), 0.0f, 1.0f)*2.0f*1.66676f; 
            }

            const float dM = pow(2.0f, 7.0f);

            vec3 CurrentDensity = vec3(dM * (PointDensity * PointDensity)) * ScatteringColor;
            vec3 S = SampleSigmaS * CurrentDensity;
            float Transmittance = exp(-(SigmaS * AirDensity) * StepSize); // exponential volumetric transmittance function
			OutputTransmittance *= Transmittance;
            vec3 IntegratedScattering = (S - S * Transmittance) / (SigmaS * AirDensity);  // S-St/Se
            // vec3 IntegratedScattering = max(exp(-PointDensity * SigmaT), exp(-PointDensity * SigmaT * 0.2f) * 0.75f) * SigmaT * (max(0.5f, SigmaS));
            VolumetricLighting += IntegratedScattering;

            // exp2() looks more correct.. somehow?
            //CurrentTransmittance *= exp2(AirDensity);
		    //TotalTransmittance *= exp2(-AirDensity);
        }

        RayPositionActual += RayDirection * StepSize;
    }

    const float IsotropicScatter = 0.25f / PI;

    /*
    vec3 ScatteredTotal  = VolumetricLighting * IsotropicScatter;
	ScatteredTotal += AmbientScatter * Ambient * IsotropicScatter;
	ScatteredTotal *= FogColor * 30.0 * (1.0 - TotalTransmittance) * StepSize;
    */

    const float ScattersM = pow(2.0f, 4.0f);
    VolumetricLighting *= ScattersM * u_Strength;
    o_Volumetrics = vec4(clamp(VolumetricLighting,0.0f,4.0f), clamp(OutputDensitySum,0.0f,0.9f/1.0f)); 

    if (false) {
        o_Volumetrics = vec4(clamp(VolumetricLighting,0.0f,5.0f), clamp(OutputTransmittance,0.0f,0.9f/1.0f)); 
    }
}


