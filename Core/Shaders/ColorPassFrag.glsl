// Light combine / Color pass.

#version 450 core

#define CLOUD_HEIGHT 70
#define PI 3.14159265359
#define THRESH 1.41414

#define clamp01(x) (clamp(x, 0.0f, 1.0f))
#define square(x) (x * x)

// Outputs 
layout (location = 0) out vec3 o_Color;
layout (location = 1) out float o_EmissivityData;

// Vertex shader inputs 
vec2 g_TexCoords;
in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;
in vec2 v_CloudSourceOcclusion;

// Uniforms
uniform sampler2D u_DiffuseTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InitialTracePositionTexture;
uniform sampler2D u_BlockIDTexture;
uniform sampler2D u_ShadowTexture;

uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlueNoiseTextures;
uniform sampler2DArray u_BlockEmissiveTextures;

uniform sampler2D u_GBufferAlbedos;
uniform sampler2D u_GBufferNormals;
uniform sampler2D u_GBufferPBR;
uniform sampler2D u_GBufferAO;




uniform samplerCube u_Skybox;
uniform sampler2D u_CloudData;
uniform sampler2D u_PreviousNormalTexture; 
uniform sampler2D u_VXAO;
uniform sampler2D u_BlueNoiseHighRes;

uniform bool u_UseDFG;
uniform bool u_RemoveTiling;
uniform bool u_CloudCatmullRomUpsampling;

uniform bool u_DEBUGDiffuseGI;
uniform bool u_DEBUGSpecGI;
uniform bool u_DEBUGShadows;


uniform sampler2D u_SSSShadowMap;
uniform bool u_SSSSS;
uniform float u_SubsurfaceScatterStrength;



uniform sampler2D u_DiffuseSHy;
uniform sampler2D u_DiffuseCoCg;

uniform float u_TextureDesatAmount;

uniform sampler2D u_SpecularIndirect;

uniform sampler2D u_HighResBL;

// Light Propogation Volume debug stuff
uniform sampler3D u_LPVLightLevel;
uniform usampler3D u_LPVColorData;
uniform int u_LPVDebugState; // Debug state 


uniform sampler2D u_DebugTexture; 


uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;
uniform vec3 u_StrongerLightDirection;
uniform float u_SunStrengthModifier;
uniform float u_MoonStrengthModifier;

uniform float u_Time;
uniform float u_GrassblockAlbedoID;
uniform float u_POMHeight;
uniform bool u_DitherPOM;

uniform mat4 u_ShadowView;
uniform mat4 u_ShadowProjection;
uniform mat4 u_ReflectionView;
uniform mat4 u_ReflectionProjection;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform vec3 u_ViewerPosition;
uniform vec2 u_Dimensions;
uniform vec2 u_CameraPlanes;

uniform bool u_CloudsEnabled;
uniform bool u_POM = false;
uniform bool u_HighQualityPOM = false;
uniform bool u_RTAO;
uniform bool u_ContactHardeningShadows;
uniform bool u_AmplifyNormalMap;
uniform bool u_VXAOCutoff;

uniform bool u_DoVXAO = true;
uniform bool u_SVGFEnabled = true;
uniform bool u_ShouldDitherUpscale = true;
uniform bool u_InferSpecularDetailSpatially = false;

uniform float u_CloudBoxSize;

uniform vec2 u_Halton;

uniform int u_LavaBlockID;
uniform sampler3D u_LavaTextures[2];

uniform samplerCube u_NebulaLowRes;
uniform float u_NebulaStrength;
uniform bool u_NebulaCelestialColor;


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

// SSBOS
layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
	int BlockSSSSSData[128];
};

layout (std430, binding = 1) buffer SSBO_BlockAverageData
{
    vec4 BlockAverageColorData[128]; // Returns the average color per block type 
};
//

// Functions :
vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 radiance_s, vec3 albedo, vec3 normal, vec3 pbr, float shadow);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
vec3 GetSmoothLPVData(vec3 UV);
vec3 GetSmoothLPVDensity(vec3 UV);
bool CompareFloatNormal(float x, float y);
vec3 GetNormalFromID(float n);
vec3 SampleNormal(sampler2D samp, vec2 txc);
vec3 GetRayDirectionAt(vec2 screenspace);
vec4 SamplePositionAt(sampler2D pos_tex, vec2 txc);
int GetBlockID(vec2 txc);
bool InScreenSpace(vec2 x);
vec2 hash2();
vec3 SampleCone(vec2 Xi, float CosThetaMax); 
vec3 XYZToRGB(in vec3 xyz);
vec3 RGBToXYZ(in vec3 rgb);
float MiePhaseFunction(float x, float g);
vec3 FresnelSchlickRoughness(vec3 Eye, vec3 norm, vec3 F0, float roughness);
vec3 FresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness);
float FastDistance(vec3 p1, vec3 p2);
vec4 BetterTexture(sampler2D samp, vec2 uv);
vec3 BasicSaturation(vec3 Color, float Adjustment);
float RayCapsuleIntersection(in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float r);
vec3 saturate(vec3 x);
float GetLuminance(vec3 color);
vec4 ClampedTexture(sampler2D tex, vec2 txc);





// Constants / Globals 


const vec3 ATMOSPHERE_SUN_COLOR = vec3(1.0f);
const vec3 ATMOSPHERE_MOON_COLOR =  vec3(0.1f, 0.1f, 1.0f);

vec3 SAMPLED_SUN_COLOR = vec3(0.0f);
vec3 SAMPLED_MOON_COLOR = vec3(0.0f);
vec3 SAMPLED_MOON_COLOR_RAW = vec3(0.0f);
vec3 SAMPLED_COLOR_MIXED = vec3(0.0f);
float HASH2SEED = 0.0f;




// Fetches atmosphere and calculates celestial spheres ->

bool GetAtmosphere(inout vec3 atmosphere_color, in vec3 in_ray_dir, float transmittance, float true_transmittance)
{
    vec3 sun_dir = (u_SunDirection); 
    vec3 moon_dir = vec3(-sun_dir.x, -sun_dir.y, sun_dir.z); 

    vec3 ray_dir = (in_ray_dir);
    
    if (true_transmittance > 0.325f) {

        if(dot(ray_dir, sun_dir) > 0.999825f)
        {
            atmosphere_color = ATMOSPHERE_SUN_COLOR; 
            o_EmissivityData = -8.0f;
           // o_BloomAlbedos = vec3(3.0f,3.0f,2.0f);
            return true;
        }
        
        if(dot(ray_dir, moon_dir) > 0.99986f)
        {
            atmosphere_color = ATMOSPHERE_MOON_COLOR;
            o_EmissivityData = -16.0f;
           // o_BloomAlbedos = vec3(0.6f,0.6f,1.0f) * 1.0f;
            return true;
        }
    }

    vec3 atmosphere = texture(u_Skybox, ray_dir).rgb;
    atmosphere_color = atmosphere;

    return false;
}



// fetch sky
vec3 GetAtmosphereAndClouds()
{
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; 

    vec3 NormalizedDir = normalize(v_RayDirection);

    vec2 SampleCoord = g_TexCoords;

	vec4 SampledCloudData;

    if (u_CloudCatmullRomUpsampling) {
        //SampledCloudData = texture_catmullrom(u_CloudData, SampleCoord).rgba;
        SampledCloudData = BetterTexture(u_CloudData, SampleCoord).rgba;
	}
    
    else {
        SampledCloudData = BetterTexture(u_CloudData, SampleCoord).rgba;
    }


    vec3 Sky = vec3(0.0f);

    bool v = GetAtmosphere(Sky, NormalizedDir, SampledCloudData.w * 20.5f, SampledCloudData.w);

    // Night sky color shifts ->
    Sky = mix(BasicSaturation(Sky, 0.8f), Sky, SunVisibility);
    Sky *= mix(2.6f, 1.0f, SunVisibility);

    if (!u_CloudsEnabled) {

        return Sky;
    }

    float transmittance = max(SampledCloudData.w, 0.01f);
    return (Sky * transmittance) + SampledCloudData.xyz;

}

vec3 FetchSky() { return GetAtmosphereAndClouds(); }


bool GetPlayerIntersect(in vec3 WorldPos, in vec3 d)
{
	float x = 0.5 - 0.3f;
	vec3 VP = u_ViewerPosition + vec3(-x, -x, +x);
    float t = RayCapsuleIntersection(WorldPos, d, VP, VP + vec3(0.0f, 1.0f, 0.0f), 0.5f);
    return t > 0.0f;
}



int RNG_SEED;

float ComputeShadow(vec3 world_pos, vec3 flat_normal, float upscaled)
{
    vec2 Txc = v_TexCoords;//ReprojectShadow(world_pos);
    float BaseShadow = upscaled;//texture(u_ShadowTexture, Txc).r;
	
    float Shadow = BaseShadow;
    float PlayerShadow  = 0.0f;
    int PlayerShadowSamples = 6;

    vec3 BiasedWorldPos = world_pos - (flat_normal * 0.4f);
    if (u_ContactHardeningShadows) {
		vec3 L = (u_StrongerLightDirection);
		vec3 T = normalize(cross(L, vec3(0.0f, 1.0f, 1.0f)));
		vec3 B = cross(T, L);
		mat3 TBN = mat3(T, B, L);
	
        for (int i = 0 ; i < PlayerShadowSamples ; i++)
        {
            vec2 Hash = hash2();
	        const float CosTheta = 0.9999905604617f; // -> changed to reduce variance. THIS IS NOT PHYSICALLY CORRECT
			vec3 ConeSample = TBN * SampleCone(Hash.xy, CosTheta); 
            PlayerShadow += float(GetPlayerIntersect(BiasedWorldPos.xyz, ConeSample.xyz));
        }

        PlayerShadow /= float(PlayerShadowSamples);
    }

    else 
    {
        PlayerShadow = float(GetPlayerIntersect(BiasedWorldPos.xyz, (u_StrongerLightDirection).xyz));
    }

    return clamp(BaseShadow + PlayerShadow, 0.0f, 1.0f);
}


bool IsInScreenSpaceBounds(in vec2 tx)
{
    if (tx.x > 0.0f && tx.y > 0.0f && tx.x < 1.0f && tx.y < 1.0f)
    {
        return true;
    }

    return false;
}


float SHToY(vec4 shY)
{
    return max(0, 3.544905f * shY.w); // get luminance (Y) from the first spherical harmonic
}

// from quake2rtx
vec3 SHToIrridiance(vec4 shY, vec2 CoCg, vec3 v)
{
    float x = dot(shY.xyz, v);
    float Y = 2.0 * (1.023326f * x + 0.886226f * shY.w);
    Y = max(Y, 0.0);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

vec3 SHToIrridiance(vec4 shY, vec2 CoCg)
{
    float Y = max(0, 3.544905f * shY.w);
    Y = max(Y, 0.0);
	CoCg *= Y * 0.282095f / (shY.w + 1e-6);
    float T = Y - CoCg.y * 0.5f;
    float G = CoCg.y + T;
    float B = T - CoCg.x * 0.5f;
    float R = B + CoCg.x;
    return max(vec3(R, G, B), vec3(0.0f));
}

// Approximates the environment brdf integration map
// Quick to sample, fixes some sampling artifacts 
vec2 KarisEnvBRDFApprox(float NdotV, float roughness)
{
	vec4 c0 = vec4(-1., -0.0275, -0.572, 0.022);
	vec4 c1 = vec4(1., 0.0425, 1.040, -0.040);
	vec4 r = roughness * c0 + c1;
	float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
	return vec2(-1.04, 1.04) * a004 + r.zw;
}

// Cook Torrance brdf : 
float ndfGGX(float cosLh, float roughness)
{
	float alpha = roughness * roughness;
	float alphaSq = alpha * alpha;
	float denom = (cosLh * cosLh) * (alphaSq - 1.0) + 1.0;
	return alphaSq / (PI * denom * denom);
}

float gaSchlickG1(float cosTheta, float k)
{
	return cosTheta / (cosTheta * (1.0 - k) + k);
}

float gaSchlickGGX(float cosLi, float cosLo, float roughness)
{
	float r = roughness + 1.0;
	float k = (r * r) / 8.0; // Epic suggests using this roughness remapping for analytic lights.
	return gaSchlickG1(cosLi, k) * gaSchlickG1(cosLo, k);
}

vec3 fresnelSchlick(vec3 F0, float cosTheta)
{
	return F0 + (vec3(1.0) - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 radiance_s, vec3 albedo, vec3 normal, vec3 pbr, float shadow)
{
    const float DIELECTRIC = 0.04f;
    const float Epsilon = 0.00001;
    float Shadow = min(shadow, 1.0f);

    vec3 Lo = normalize(u_ViewerPosition - world_pos);

	vec3 N = normal;
	float cosLo = max(0.0, dot(N, Lo));
	vec3 Lr = 2.0 * cosLo * N - Lo;
	vec3 F0 = mix(vec3(DIELECTRIC), albedo, pbr.g);

    vec3 Li = light_dir;
	vec3 Lradiance = radiance;

	vec3 Lh = normalize(Li + Lo);

	float cosLi = max(0.0, dot(N, Li));
	float cosLh = max(0.0, dot(N, Lh));

	//vec3 F  = fresnelSchlick(F0, max(0.0, dot(Lh, Lo)));
	vec3 F = FresnelSchlickRoughness(Lo, normal.xyz, vec3(F0), pbr.r); // use fresnel schlick roughness, approximates better.
	float D = ndfGGX(cosLh, pbr.r);
	float G = gaSchlickGGX(cosLi, cosLo, pbr.r);

	vec3 kd = mix(vec3(1.0) - F, vec3(0.0), pbr.g);
	vec3 diffuseBRDF = kd * albedo;

	vec3 specularBRDF = (F * D * G) / max(Epsilon, 4.0 * cosLi * cosLo); specularBRDF = clamp(specularBRDF, 0.0f, 2.0f);
	vec3 Result = (diffuseBRDF * Lradiance * cosLi) + (specularBRDF * radiance_s * cosLi);
    return clamp(Result, 0.0f, 2.5) * clamp((1.0f - Shadow), 0.0f, 1.0f);
}

float G_Smith_over_NdotV(float roughness, float NdotV, float NdotL)
{
    float alpha = square(roughness);
    float g1 = NdotV * sqrt(square(alpha) + (1.0 - square(alpha)) * square(NdotL));
    float g2 = NdotL * sqrt(square(alpha) + (1.0 - square(alpha)) * square(NdotV));
    return 2.0 *  NdotL / (g1 + g2);
}

// ggx specular
float SpecularGGX(vec3 V, vec3 L, vec3 N, float roughness, float NoH_offset)
{
    vec3 H = normalize(L - V);
    float NoL = max(0, dot(N, L));
    float VoH = max(0, -dot(V, H));
    float NoV = max(0, -dot(N, V));
    float NoH = clamp(dot(N, H) + NoH_offset, 0, 1);

    if (NoL > 0)
    {
        float G = G_Smith_over_NdotV(roughness, NoV, NoL);
        float alpha = square(max(roughness, 0.02));
        float D = square(alpha) / (PI * square(square(NoH) * square(alpha) + (1 - square(NoH))));
        return D * G / 4;
    }

    return 0;
}

// Custom very non physically based subsurface scattering 
vec3 IntegrateSubsurfaceScatter(vec3 V, vec3 P, vec3 N, vec3 Albedo, float Shadow, float SunVisibility) 
{
    float VDotL = dot(V, u_StrongerLightDirection);
    float MiePhase = MiePhaseFunction(VDotL, 0.8f);
    MiePhase = pow(MiePhase,(0.9f/1.0f)*PI*(0.9f/1.0f));
    float ScatteringAmount = (1.0f - exp(-pow(MiePhase * PI * 2.0f * PI, 1.0f)));
    ScatteringAmount = clamp(ScatteringAmount * PI * 1.1f * u_SubsurfaceScatterStrength * mix(0.5f, 1.0f, SunVisibility), 0.0001f, 7.0f);
    float VisibilityTerm = (1.0f - Shadow);
    vec3 Scattering = vec3(ScatteringAmount * mix(SAMPLED_MOON_COLOR,SAMPLED_SUN_COLOR/3.5f,SunVisibility) * VisibilityTerm) * Albedo;
    return Scattering;
}

// Derives approximate specular from diffuse spherical harmonic 
// credits : q2rtx
vec3 DeriveApproxSpecular(vec4 SHy, vec3 IndirectDiffuse, vec3 Eye, vec3 Normal, float rough) 
{  
	float Roughness = clamp(rough, 0.1f, 0.45f);
	vec3 IncomingDir = SHy.xyz / SHy.w * (0.282095f / 0.488603f);
	vec3 RawSpecularDir = reflect(Eye, Normal); 
	float IncomingLen = length(IncomingDir);
	float Directionality = IncomingLen;
    float Scale = 1.0f;

	if(Directionality >= 1.0) 
    {
		IncomingDir /= IncomingLen; 
	}

	else
	{
		IncomingDir = mix(RawSpecularDir, IncomingDir / (IncomingLen + 0.00001f), vec3(Directionality));
		Scale = pow(Roughness + 1.0f, 3.0f);
	}

	float SpecularGGXIntegrated = SpecularGGX(Eye, IncomingDir, Normal, max(Roughness, 0.01f), 0.0f);
	vec3 Integrated = pow(SpecularGGXIntegrated, 1.2f) * IndirectDiffuse * 18.0f * Scale;
	if (isnan(Integrated.x)||isinf(Integrated.x)||isnan(Integrated.y)||isinf(Integrated.y)||isnan(Integrated.z)||isinf(Integrated.z)) { Integrated = vec3(0.0f); }
	return max(Integrated, 0.00001f); 

}


bool IsImmediateNeighbour(ivec2 x) {
    if (x == ivec2(-1, 0) || x == ivec2(1, 0) || x == ivec2(0, -1) || x == ivec2(0, 1) || x == ivec2(0, 0)) {
        return true;
    }

    return false;
}


// -- Spatial upscaler based on atrous wavelet filter --

void SpatiallyUpscaleBuffers(vec3 BaseNormal, float BaseLinearDepth, out vec4 SH, out vec2 CoCg, out vec4 SpecularIndirect, out float ShadowSample, bool fuckingsmooth, out float ao)
{
    const bool DoSpatialUpscaling = true;

    if (! DoSpatialUpscaling) {
        vec2 SampleCoord = g_TexCoords;
        SH = texture(u_DiffuseSHy, SampleCoord).xyzw;
		CoCg += texture(u_DiffuseCoCg, SampleCoord).xy;
        SpecularIndirect += texture(u_SpecularIndirect, SampleCoord).xyzw;
        ShadowSample += texture(u_ShadowTexture, SampleCoord).x;
        ao = texture(u_VXAO, g_TexCoords).x;
		
        return;
    }

	vec4 TotalSH = vec4(0.0f);
	vec2 TotalCoCg = vec2(0.0f);
    vec4 TotalSpecColor = vec4(0.0f);
    float TotalShadow = 0.0f;
	float TotalWeight = 0.0f;
	float TotalShadowWeight = 0.0f;
	float TotalSpecWeight = 0.0f;
    ao = 0.0f;

    bool FilterShadows = BaseLinearDepth > 128.0f;

    if (!FilterShadows) {
        TotalShadow = BetterTexture(u_ShadowTexture, g_TexCoords).x;
    }

	vec2 TexelSize = 1.0f / textureSize(u_DiffuseSHy, 0);
    const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

	for (int x = -1 ; x <= 1 ; x++) {

        for (int y = -1 ; y <= 1 ; y++) {
            
			// Calculate sample coordinate ->
            vec2 SampleCoord = g_TexCoords + (vec2(x,y)) * TexelSize;
			
			// Depth weight ->
		    float LinearDepthAt = (1.0f/texture(u_InitialTracePositionTexture, SampleCoord).x);
            float ExpDepth = exp(-(abs(LinearDepthAt - BaseLinearDepth) * 2.0f));
            float DepthWeight = pow(ExpDepth, 5.25f);

			// Normal weight ->
            vec3 NormalAt = SampleNormal(u_NormalTexture, SampleCoord.xy).xyz;
		    float NormalWeight = pow(max(dot(NormalAt, BaseNormal),0.000001f), 16.0f);
            NormalWeight = max(NormalWeight, 0.0001f);

			// Kernel weight ->
            float KernelWeight = AtrousWeights[abs(x)] * AtrousWeights[abs(y)];

			// Combine weights ->
            float Weight = clamp(NormalWeight * DepthWeight * KernelWeight, 0.0f, 1.0f);
           
			// Specular weight ->
		    float SpecWeight = Weight;//clamp(NormalWeight * NormalWeight * pow(ExpDepth, mix(4.25f, 12.0f, BaseLinearDepth < 32.0f)) * KernelWeight, 0.0f, 1.0f);
		    
			// To avoid division by 0
			Weight = max(Weight, 0.00000001f);
		    SpecWeight = max(SpecWeight, 0.00000001);

			// Sum ->
			
			// Diffuse 
		    TotalSH += texture(u_DiffuseSHy, SampleCoord).xyzw * Weight;
		    TotalCoCg += texture(u_DiffuseCoCg, SampleCoord).xy * Weight;
			
			// Spec 
            TotalSpecColor += texture(u_SpecularIndirect, SampleCoord).xyzw * SpecWeight;
			
			// Ambient occlusion 
            ao += texture(u_VXAO, SampleCoord).x * Weight;
			
			// Weights 
		    TotalWeight += Weight;
            TotalSpecWeight += SpecWeight;

			// Filter shadows if current offset is an immediate neighbour
            if (FilterShadows) {
				bool SampleDirect = IsImmediateNeighbour(ivec2(x,y));
				
                if (SampleDirect) {
                    float ShadowWeight = KernelWeight * clamp(pow(ExpDepth, 32.0f), 0.00001f, 1.0f) * (NormalWeight * NormalWeight);
		            TotalShadow += texture(u_ShadowTexture, SampleCoord).x * ShadowWeight;
                    TotalShadowWeight += ShadowWeight;
                }
            }
        }
	}

    TotalWeight = max(TotalWeight, 0.0000001f);

	
    SH = TotalSH / TotalWeight;
    CoCg = TotalCoCg / TotalWeight;
    TotalSpecColor = TotalSpecColor / TotalSpecWeight;
    ao = ao / TotalWeight;

    if (FilterShadows) {
        TotalShadow = TotalShadow / TotalShadowWeight;
    }


	// We use low precision buffers for sun shadows, this is to avoid precision artifacts (banding) ->
	TotalShadow  += 0.05f;
	
	TotalShadow = clamp(TotalShadow, 0.00000000000001f, 1.0f);

    SpecularIndirect = TotalSpecColor;
    ShadowSample = TotalShadow;
	
}


// Converts a temperature to a color
vec3 TemperatureToRGB(float temperatureInKelvins);

vec3 SetColorLuminance(vec3 Color, float L) {
    Color = RGBToXYZ(Color);
    Color.y = L;
    Color = XYZToRGB(Color);
    return Color;
}

vec3 SetColorLuminance(vec3 Color, vec3 Color2) {
    Color = RGBToXYZ(Color);
    Color.y = RGBToXYZ(Color2).y;
    Color = XYZToRGB(Color);
    return Color;
}


vec3 SampleSunColor()
{
    const vec3 TemperatureModifier = TemperatureToRGB(5778.0f);
    vec3 SunTransmittance = texture(u_Skybox, u_SunDirection.xyz).xyz;
    vec3 SunColor = SunTransmittance;
    SunColor *= TemperatureModifier;
    return SunColor * PI * 2.2f * u_SunStrengthModifier;
}

vec3 SampleMoonColor(float sv, bool x) {
    vec3 MoonTransmittance = texture(u_Skybox, u_MoonDirection).xyz;
    vec3 MoonColor = MoonTransmittance;

    if (u_NebulaCelestialColor)
    {
        float TimeOfDayTransition = sv * sv;
        vec3 Nebula = BasicSaturation(texture(u_NebulaLowRes, u_MoonDirection).xyz, 1.25f) * 0.5f * u_NebulaStrength * TimeOfDayTransition;
        float MoonTransmittanceY = RGBToXYZ(MoonTransmittance).y;

        // Clamp luminance ->
        MoonColor = MoonTransmittance + Nebula;
        MoonColor = RGBToXYZ(MoonColor);
        MoonColor.y = mix(MoonColor.y, MoonTransmittanceY, 0.8f);
        MoonColor = XYZToRGB(MoonColor);

        if (!x) {
            MoonColor = mix(MoonColor, MoonTransmittance, 0.9f);
        }
    }

    MoonColor = MoonColor * PI * u_MoonStrengthModifier;
    float L = GetLuminance(MoonColor);
    //MoonColor = mix(MoonColor, vec3(L), 0.05f); 
    return MoonColor * 0.5f * vec3(0.9f, 0.9f, 1.0f);
}

// Converts a value from one range to another
float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

vec3 LavaTextureDistortion(vec3 UV) {
    vec2 UVxy = UV.xy;
    float time = u_Time;
    UV.x += sin(time * 0.25f);
    UV.y += pow(cos(time * 0.15f),2.);
    UV.x += cos(UV.x*10.0f + time)*0.3f;
    UV.y += sin(UV.y*5.0f + UV.x*4.0f + time*1.3f)*0.4f;
    UV.xy = mix(UV.xy,UVxy.xy,0.91f);
    return UV;
}

float DelinearizeDepth(float d) 
{
	float near = u_CameraPlanes.x;
	float far = u_CameraPlanes.y;
    d = d < 0.0f ? d : d + near;
	d = d < 0.0f ? far - 0.01f : d;
	float Delinearized = -((near + far) * d - (2.0f * near)) / ((near - far) * d);
	return Delinearized;
}


float ComputeSpecularAO(float NoV, float ao, float roughness) 
{
    return clamp(pow(NoV + ao, exp2(-16.0f * roughness - 1.0f)) - 1.0f + ao, 0.0f, 1.0f);
}

void main()
{
    g_TexCoords = v_TexCoords;
    //g_TexCoords += u_Halton/u_Dimensions;
    o_EmissivityData.x = (-1.0f);


	RNG_SEED = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(100.0f * fract(u_Time));

    // Xorshift!
	RNG_SEED ^= RNG_SEED << 13;
    RNG_SEED ^= RNG_SEED >> 17;
    RNG_SEED ^= RNG_SEED << 5;
    HASH2SEED = (g_TexCoords.x * g_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(u_Time) * 100.0f;

    vec4 WorldPosition = SamplePositionAt(u_InitialTracePositionTexture, g_TexCoords);

    vec3 SampledNormals = SampleNormal(u_NormalTexture, g_TexCoords).rgb;
    vec3 AtmosphereAt = vec3(0.0f);
    o_Color = vec3(1.0f);

    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
    SAMPLED_SUN_COLOR = SampleSunColor();
    SAMPLED_MOON_COLOR = SampleMoonColor(SunVisibility, true);
    SAMPLED_MOON_COLOR_RAW = SampleMoonColor(SunVisibility, false);
    SAMPLED_COLOR_MIXED = mix(SAMPLED_SUN_COLOR, SAMPLED_MOON_COLOR, SunVisibility);
    int BaseBlockID = GetBlockID(g_TexCoords);

    if (WorldPosition.w > 0.0f)
    {
       

        float Bias = 0.0035f;
        bool InBiasedSS =  (g_TexCoords.x > 0.0 + Bias && g_TexCoords.x < 1.0 - Bias 
		 && g_TexCoords.y > 0.0 + Bias && g_TexCoords.y < 1.0f - Bias);

        //if (!IsAtEdge(g_TexCoords))
        if (true)
        {
            bool DoSSS = BlockSSSSSData[BaseBlockID] > 0 && u_SSSSS;






            bool IsLiquid =  (u_LavaBlockID == BaseBlockID) ;
            bool IsLava =  (u_LavaBlockID == BaseBlockID) ;

            vec3 AlbedoColor = texture(u_GBufferAlbedos, g_TexCoords).xyz;
                               
            vec3 NormalMapped = texture(u_GBufferNormals, g_TexCoords).xyz;

            vec4 PBRMap = texture(u_GBufferPBR, g_TexCoords).xyzw;
            float Emissivity = PBRMap.w;
            PBRMap.w = texture(u_GBufferAO, g_TexCoords).x;

            AlbedoColor = BasicSaturation(AlbedoColor, 1.0f - u_TextureDesatAmount);

            if (PBRMap.y >= 0.1f - 0.01f){
                AlbedoColor = BasicSaturation(AlbedoColor, 0.9f);
            }

            vec3 NonAmplifiedNormal = NormalMapped;

            if (u_AmplifyNormalMap) {
                NormalMapped.x *= 1.64f;
                NormalMapped.z *= 1.85f;
                NormalMapped += 1e-4f;
                NormalMapped = normalize(NormalMapped);
            }
            
            Emissivity *= mix(1.,1.5,float(IsLava));

            vec4 SHy;
            vec2 ShCoCg;
            vec4 UpscaledSpecularIndirect;
			float UpscaledShadow=0.0f;
            float UpscaledAO = 0.0f;
		
            SpatiallyUpscaleBuffers(SampledNormals.xyz, WorldPosition.w, SHy, ShCoCg, UpscaledSpecularIndirect,UpscaledShadow, PBRMap.x <= 0.1, UpscaledAO);
			UpscaledShadow = clamp(UpscaledShadow,0.0f,1.0f);

			float RayTracedShadow = ComputeShadow(WorldPosition.xyz, SampledNormals.xyz,UpscaledShadow);
			RayTracedShadow = clamp(RayTracedShadow,0.,1.);


            vec3 IndirectN = NormalMapped.xyz;
            vec3 SampledIndirectDiffuse = vec3(0.0f);

            if (true) {

                //vec3 MostProminentDirectionDiffuse = SHy.xyz / SHy.w * (0.282095f / 0.488603f);
                SampledIndirectDiffuse = SHToIrridiance(SHy, ShCoCg, IndirectN);

           
            }


            bool CorrectedSpecular = false;

            
            // VXAO from indirect diffuse trace : 

            bool do_vxao = u_DoVXAO && u_SVGFEnabled;
			
			const int VXGI_CUTOFF = 196;
			
            if (do_vxao)
            {
                float bias = 0.00125f;
                if (g_TexCoords.x > bias && g_TexCoords.x < 1.0f - bias &&
                    g_TexCoords.y > bias && g_TexCoords.y < 1.0f - bias)
                {
                    float ao = UpscaledAO;//BetterTexture(u_VXAO, g_TexCoords).x;

                    //float fade = u_VXAOCutoff ? (1.0f - exp(-WorldPosition.w * 0.00125f)) : 0.0f;
                    float fade = u_VXAOCutoff ? (WorldPosition.w < VXGI_CUTOFF ? 0.0f : 1.0f) : 0.0f;

                    ao = mix(ao, 1.0f, fade);
					
                    SampledIndirectDiffuse.xyz *= vec3(clamp(pow(ao, 2.2f), 0.05, 1.0f));

                    //o_Color = vec3(pow(ao, 4.2f)*0.7f);
                    //return;
                }
            }

            vec3 LightAmbience = (vec3(120.0f, 172.0f, 255.0f) / 255.0f) * 1.01f;
            vec3 Ambient = (AlbedoColor * LightAmbience) * 0.09f;
            float SampledAO = pow(PBRMap.w, 1.25f);

            vec3 SunColor = SAMPLED_SUN_COLOR;
            vec3 MoonColor = SAMPLED_MOON_COLOR;

            vec3 SunDirectLighting = CalculateDirectionalLight(WorldPosition.xyz, (u_SunDirection), SunColor, SunColor, AlbedoColor, NormalMapped, PBRMap.xyz, RayTracedShadow);
            vec3 MoonDirectLighting = CalculateDirectionalLight(WorldPosition.xyz, (u_MoonDirection), MoonColor, MoonColor, AlbedoColor, NormalMapped, PBRMap.xyz, RayTracedShadow);
            vec3 DirectLighting = mix(SunDirectLighting, MoonDirectLighting, SunVisibility * vec3(1.0f));
            
            DirectLighting = (float(!(Emissivity > 0.05f)) * DirectLighting);
			
			DirectLighting = max(DirectLighting, 0.000001f);
			
            float Roughness = PBRMap.r;
			
            vec3 DiffuseIndirect = u_LPVDebugState == 3 ? ((GetSmoothLPVDensity(WorldPosition.xyz+SampledNormals.xyz) * AlbedoColor)+(texture(u_VXAO,v_TexCoords).y*BasicSaturation(texture(u_Skybox,vec3(0.0f,1.0f,0.0f)).xyz,0.75f)*3.*AlbedoColor)) : 
                                  (u_LPVDebugState == 4 ? ((GetSmoothLPVData(WorldPosition.xyz+SampledNormals.xyz) * AlbedoColor * 5.0f)+(texture(u_VXAO,v_TexCoords).y*BasicSaturation(texture(u_Skybox,vec3(0.0f,1.0f,0.0f)).xyz,0.75f)*3.*AlbedoColor)) : (SampledIndirectDiffuse.xyz * AlbedoColor));
           
            vec3 SpecularIndirect = vec3(0.0f);
            
            vec3 I = normalize(WorldPosition.xyz - u_ViewerPosition);
            vec3 R = reflect(I, NormalMapped).xyz;

            if (InBiasedSS) {
                float Attenuation = clamp(pow(max(dot(SampledNormals.xyz, NonAmplifiedNormal), 0.00000001f),2.0f)*1.5f, 0.1f, 1.0f);
                SpecularIndirect = UpscaledSpecularIndirect.xyz * 1.5f * Attenuation;
            }

			SpecularIndirect = max(vec3(0.0001f), SpecularIndirect);
            


            vec3 Lo = normalize(u_ViewerPosition - WorldPosition.xyz); // Outgoing direction 
            vec3 F0 = mix(vec3(0.04), AlbedoColor, PBRMap.g); // Fresnel at 0 degrees.

            
            // Final combine : 
            // Use fresnel to get the amount of diffuse and reflected light 
            DiffuseIndirect *= clamp(SampledAO, 0.2f, 1.01f);

            bool dfg = u_UseDFG;

            vec3 SubsurfaceScatter = vec3(0.);

            if (DoSSS) {
				const float FakeRefractiveIdx = 1.0f/1.333333f;
				vec3 FakeRefracted = (refract(-Lo,NormalMapped,FakeRefractiveIdx));
				float SSSShadowMapFetch = texture(u_SSSShadowMap, v_TexCoords).x;
                SubsurfaceScatter = IntegrateSubsurfaceScatter(FakeRefracted, WorldPosition.xyz, SampledNormals.xyz, AlbedoColor, SSSShadowMapFetch, 1.-SunVisibility);
            }

            if (!dfg) {
                vec3 SpecularFactor = FresnelSchlickRoughness(Lo, NormalMapped.xyz, vec3(F0), Roughness); 
                
                if (PBRMap.y >= 0.1f) {
                    SpecularIndirect *= 1.1f;
                }

                o_Color = ((DirectLighting + SubsurfaceScatter) + ((1.0f - SpecularFactor) * DiffuseIndirect) + 
                      (SpecularFactor * SpecularIndirect));
            }

            else {

                
				
				// Account for Specular AO ->
				float NoV = clamp(dot(NormalMapped, Lo), 0.0000000001f, 1.0f);
				float SpecularAO = ComputeSpecularAO(NoV, PBRMap.w, Roughness);
				SpecularAO = pow(SpecularAO, 1.5f);
				
				// Light leak handling -> (from filament BRDF) 
				float HorizonAmount = min(1.0f + dot(R, NormalMapped), 1.0f);
				HorizonAmount *= HorizonAmount * HorizonAmount;
				

				// Not physically accurate
                // Done for stylization purposes 
                // Metals have their reflections 60% brighter and have their albedos 10% more desaturated
                
				
				if (PBRMap.y >= 0.1f) {
                    SpecularIndirect *= 1.6f;
                }
				
				SpecularIndirect *= SpecularAO;
				SpecularIndirect *= HorizonAmount;

                vec3 FresnelTerm = FresnelSchlickRoughness(Lo, NormalMapped.xyz, vec3(F0), Roughness); 
                FresnelTerm = clamp(FresnelTerm, 0.0f, 1.0f);
                vec3 kS = FresnelTerm;
                vec3 kD = 1.0f - kS;
                kD *= 1.0f - PBRMap.y;
				
				// Environment brdf (Karis)
                vec2 EnvironmentBRDFSampleLocation = vec2(max(dot(-I, NormalMapped.xyz), 0.000001f), PBRMap.x);
                EnvironmentBRDFSampleLocation = clamp(EnvironmentBRDFSampleLocation, 0.0f, 1.0f);
                vec2 EnvironmentBRDF = KarisEnvBRDFApprox(EnvironmentBRDFSampleLocation.x, EnvironmentBRDFSampleLocation.y);
				
				// Combine using environment brdf and fresnel term
                vec3 IndirectSpecularFinal = SpecularIndirect * (FresnelTerm * EnvironmentBRDF.x + EnvironmentBRDF.y);
                vec3 IndirectLighting = (kD * DiffuseIndirect) + IndirectSpecularFinal + SubsurfaceScatter;
                o_Color = DirectLighting + IndirectLighting;
            }


			// Debug ->
				if (u_LPVDebugState == 1) {
					o_Color = GetSmoothLPVDensity(WorldPosition.xyz+SampledNormals.xyz);
				}
	
				if (u_LPVDebugState == 2) {
					o_Color = GetSmoothLPVData(WorldPosition.xyz+SampledNormals.xyz);
				}
	
				if (u_DEBUGDiffuseGI) {
					o_Color = 1.0f - exp(-SampledIndirectDiffuse);
				}
	
				if (u_DEBUGSpecGI) {
					o_Color = 1.0f - exp(-SpecularIndirect);
				}
	
				if (u_DEBUGShadows) {
					o_Color = vec3(1.0f - RayTracedShadow);
				// o_Color += bayer128(gl_FragCoord.xy) / 128.0f;
				}
			//
	
	
	
            
            o_EmissivityData = Emissivity;

             
            // Flicker
            if (IsLava) {
                o_EmissivityData *= clamp(pow(sin(u_Time * 5.5f) * 0.5f + 0.5f, 1.0f/3.0f), 0.7f, 1.0f);
                o_EmissivityData *= 1.325f;
            }

            // Approximate specular ->
            const bool DebugApproximateSpecular = false;

            if (DebugApproximateSpecular) {
                o_Color = DeriveApproxSpecular(SHy, SampledIndirectDiffuse, I, NormalMapped.xyz, Roughness);
            }


            //o_Color = DiffuseIndirect;

          // o_Color = vec3(pow(PBRMap.w, 10.0f));

        }

        else 
        {   
            vec3 CloudAndSky = GetAtmosphereAndClouds();
			o_Color = (CloudAndSky);
			//o_EmissivityData.x = (-1.0f);
        }
    }

    else 
    {   
        vec3 CloudAndSky = GetAtmosphereAndClouds();
        o_Color = (CloudAndSky);
        //o_EmissivityData.x = (-1.0f);
    }

    const bool DEBUG_SHIT = false;

    if (DEBUG_SHIT) {
        //vec4 DebugData = texture(u_DebugTexture, ProjectDirection(v_RayDirection)).xyzw;
        vec4 DebugData = texture(u_DebugTexture, g_TexCoords).xyzw;
        //o_Color = max(vec3(DebugData.x / 20.0f), 0.0f);
        o_Color = max(vec3(DebugData.xyz),0.00001f);
    }




	o_Color = max(o_Color, vec3(0.000001f));
}





// Utility -> -> ->

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

// Samplers :
vec3 SampleNormal(sampler2D samp, vec2 txc) { 
	return GetNormalFromID(texture(samp, txc).x);
}



vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SamplePositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = 1.0f/texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

int GetBlockID(vec2 txc)
{
	float id = texture(u_BlockIDTexture, txc).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

bool InScreenSpace(vec2 x)
{
    return x.x < 1.0f && x.x > 0.0f && x.y < 1.0f && x.y > 0.0f;
}

// RNG 
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}


vec3 SampleCone(vec2 Xi, float CosThetaMax) 
{
    float CosTheta = (1.0 - Xi.x) + Xi.x * CosThetaMax;
    float SinTheta = sqrt(1.0 - CosTheta * CosTheta);
    float phi = Xi.y * PI * 2.0;
    vec3 L;
    L.x = SinTheta * cos(phi);
    L.y = SinTheta * sin(phi);
    L.z = CosTheta;
    return L;
}

// Utility ->
const mat3x3 xyzToRGBMatrix = mat3(
    3.1338561, -1.6168667, -0.4906146,
    -0.9787684,  1.9161415,  0.0334540,
    0.0719453, -0.2289914,  1.4052427
);

const mat3x3 rgbToXYZMatrix = mat3(
    vec3(0.5149, 0.3244, 0.1607),
    vec3(0.3654, 0.6704, 0.0642),
    vec3(0.0248, 0.1248, 0.8504)
);

vec3 XYZToRGB(in vec3 xyz) 
{
    float r = dot(xyz, xyzToRGBMatrix[0]);
    float g = dot(xyz, xyzToRGBMatrix[1]);
    float b = dot(xyz, xyzToRGBMatrix[2]);
    return vec3(r, g, b);
}

vec3 RGBToXYZ(in vec3 rgb) 
{
    float x = dot(rgb, rgbToXYZMatrix[0]);
    float y = dot(rgb, rgbToXYZMatrix[1]);
    float z = dot(rgb, rgbToXYZMatrix[2]);
    return vec3(x, y, z);
}


float MiePhaseFunction(float x, float g){
    float gg = g * g;
    return (gg * -0.25 /3.14 + 0.25 /3.14) * pow(-2.0 * (g * x) + (gg + 1.0), -1.5);
}

//vec3 FresnelSchlickRoughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
//{
//    float cosTheta = max(dot(norm, Eye), 0.0);
//    const float amplifier = 2.4f;
//    return F0 + (max(vec3(pow(1.0f - roughness, amplifier)), F0) - F0) * pow(clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
//}

vec3 FresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness)
{
    return F0 + (max(vec3(1.0-roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 FresnelSchlickRoughness(vec3 Eye, vec3 norm, vec3 F0, float roughness)
{
    float cosTheta = clamp(dot(Eye, norm), 0.00001f, 1.0f);
    return F0 + (max(vec3(1.0-roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

float FastDistance(vec3 p1, vec3 p2)
{
    return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

vec4 BetterTexture(sampler2D samp, vec2 uv)
{
    vec2 textureResolution = textureSize(samp, 0).xy;
    uv = uv * textureResolution + 0.5;
    vec2 iuv = floor(uv);
    vec2 fuv = fract(uv);
    uv = iuv + fuv * fuv * (3.0 - 2.0 * fuv);
    uv = (uv - 0.5) / textureResolution;
    return texture(samp, uv).xyzw;
}

vec3 BasicSaturation(vec3 Color, float Adjustment)
{
    const vec3 LuminosityCoefficients = vec3(0.2125f, 0.7154f, 0.0721f);
    vec3 Luminosity = vec3(dot(Color, LuminosityCoefficients));
    return mix(Luminosity, Color, Adjustment);
}

float RayCapsuleIntersection(in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float r)
{
    vec3 ba = pb - pa;
    vec3 oa = ro - pa;

    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoa = dot(ba, oa);
    float rdoa = dot(rd, oa);
    float oaoa = dot(oa, oa);

    float a = baba - bard * bard;
    float b = baba * rdoa - baoa * bard;
    float c = baba * oaoa - baoa * baoa - r * r * baba;
    float h = b * b - a * c;

    if (h >= 0.0)
    {
        float t = (-b - sqrt(h)) / a;
        float y = baoa + t * bard;
        if (y > 0.0 && y < baba) return t;
        vec3 oc = (y <= 0.0) ? oa : ro - pb;
        b = dot(rd, oc);
        c = dot(oc, oc) - r * r;
        h = b * b - c;
        if (h > 0.0) return -b - sqrt(h);
    }
    return -1.0;
}

vec3 saturate(vec3 x)
{
    return clamp(x, 0.0f, 1.0f);
}

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

vec4 ClampedTexture(sampler2D tex, vec2 txc)
{
    return texture(tex, clamp(txc, 0.0f, 0.9999999f));
}




// LPV stuff

vec3 SampleLPVColor(vec3 UV) {
    uint BlockID = texture(u_LPVColorData, UV).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 SampleLPVColor(vec3 UV, float D) {
    uint BlockID = texture(u_LPVColorData, UV+D*0.5f*(1.0f/vec3(384.0f,128.0f,384.0f))).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);
}   

vec3 SampleLPVColorTexel(ivec3 Texel, int LOD) {
    uint BlockID = texelFetch(u_LPVColorData, Texel, LOD).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 SampleLPVColorTexel(ivec3 Texel) {
    uint BlockID = texelFetch(u_LPVColorData, Texel, 0).x;
    return vec3(BlockAverageColorData[clamp(BlockID,0,128)]);

}   

vec3 InterpolateLPVColorData(vec3 uv)
{
    vec3 res = vec3(384.0f, 128.0f, 384.0f);
    vec3 q = fract(uv * res);
    ivec3 t = ivec3(uv * res);
    ivec3 e = ivec3(-1, 0, 1);
    vec3 q0 = (q+1.0)/2.0;
    vec3 q1 = q/2.0;	
    vec3 s000 = SampleLPVColorTexel(t + e.xxx, 0);
    vec3 s001 = SampleLPVColorTexel(t + e.xxy, 0);
    vec3 s002 = SampleLPVColorTexel(t + e.xxz, 0);
    vec3 s012 = SampleLPVColorTexel(t + e.xyz, 0);
    vec3 s011 = SampleLPVColorTexel(t + e.xyy, 0);
    vec3 s010 = SampleLPVColorTexel(t + e.xyx, 0);
    vec3 s020 = SampleLPVColorTexel(t + e.xzx, 0);
    vec3 s021 = SampleLPVColorTexel(t + e.xzy, 0);
    vec3 s022 = SampleLPVColorTexel(t + e.xzz, 0);
    vec3 y00 = mix(mix(s000, s001, q0.z), mix(s001, s002, q1.z), q.z);
    vec3 y01 = mix(mix(s010, s011, q0.z), mix(s011, s012, q1.z), q.z);
    vec3 y02 = mix(mix(s020, s021, q0.z), mix(s021, s022, q1.z), q.z);
	vec3 x0 = mix(mix(y00, y01, q0.y), mix(y01, y02, q1.y), q.y);
    vec3 s122 = SampleLPVColorTexel(t + e.yzz, 0);
    vec3 s121 = SampleLPVColorTexel(t + e.yzy, 0);
    vec3 s120 = SampleLPVColorTexel(t + e.yzx, 0);
    vec3 s110 = SampleLPVColorTexel(t + e.yyx, 0);
    vec3 s111 = SampleLPVColorTexel(t + e.yyy, 0);
    vec3 s112 = SampleLPVColorTexel(t + e.yyz, 0);
    vec3 s102 = SampleLPVColorTexel(t + e.yxz, 0);
    vec3 s101 = SampleLPVColorTexel(t + e.yxy, 0);
    vec3 s100 = SampleLPVColorTexel(t + e.yxx, 0);
    vec3 y10 = mix(mix(s100, s101, q0.z), mix(s101, s102, q1.z), q.z);
    vec3 y11 = mix(mix(s110, s111, q0.z), mix(s111, s112, q1.z), q.z);
    vec3 y12 = mix(mix(s120, s121, q0.z), mix(s121, s122, q1.z), q.z);
    vec3 x1 = mix(mix(y10, y11, q0.y), mix(y11, y12, q1.y), q.y);
    vec3 s200 = SampleLPVColorTexel(t + e.zxx, 0);
    vec3 s201 = SampleLPVColorTexel(t + e.zxy, 0);
    vec3 s202 = SampleLPVColorTexel(t + e.zxz, 0);
    vec3 s212 = SampleLPVColorTexel(t + e.zyz, 0);
    vec3 s211 = SampleLPVColorTexel(t + e.zyy, 0);
    vec3 s210 = SampleLPVColorTexel(t + e.zyx, 0);
    vec3 s220 = SampleLPVColorTexel(t + e.zzx, 0);
    vec3 s221 = SampleLPVColorTexel(t + e.zzy, 0);
    vec3 s222 = SampleLPVColorTexel(t + e.zzz, 0);
    vec3 y20 = mix(mix(s200, s201, q0.z), mix(s201, s202, q1.z), q.z);
    vec3 y21 = mix(mix(s210, s211, q0.z), mix(s211, s212, q1.z), q.z);
    vec3 y22 = mix(mix(s220, s221, q0.z), mix(s221, s222, q1.z), q.z);
    vec3 x2 = mix(mix(y20, y21, q0.y), mix(y21, y22, q1.y), q.y);
    return mix(mix(x0, x1, q0.x), mix(x1, x2, q1.x), q.x);
}

vec3 InterpolateLPVColorDithered(vec3 UV) // Very few samples with fantastic results.
{ 
	const bool TemporalIntegration = true;
	vec2 OffsettedTxc = g_TexCoords + (vec2(fract(u_Time)*6., fract(u_Time)*2.)/max(vec2(0.0001f),u_Dimensions))*float(TemporalIntegration);
    vec3 Dither = texture(u_BlueNoiseHighRes, (OffsettedTxc * (u_Dimensions / vec2(textureSize(u_BlueNoiseHighRes,0).xy)))).xyz;
    const vec3 Resolution = vec3(384.0f, 128.0f, 384.0f);
    Dither /= Resolution;

    vec3 FractTexel = fract(UV * Resolution);
    vec3 LinearOffset = (FractTexel * (FractTexel - 1.0f) + 0.5f) / Resolution;
    vec3 W0 = UV - LinearOffset;
    vec3 W1 = UV + LinearOffset;
    const float DitherWeights[4] = float[4](1.0f, 1.0f, 1.0f, 1.0f);
    const float GlobalDitherNoiseWeight = 2.0f;

    vec3 Interpolated = SampleLPVColor(vec3(W0.x, W0.y, W0.z) + Dither * DitherWeights[0] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W0.y, W0.z) - Dither * DitherWeights[0] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W1.y, W0.z) + Dither * DitherWeights[1] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W0.x, W1.y, W0.z) - Dither * DitherWeights[1] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W0.x, W1.y, W1.z) + Dither * DitherWeights[2] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W1.y, W1.z) - Dither * DitherWeights[2] * GlobalDitherNoiseWeight)
    	   + SampleLPVColor(vec3(W1.x, W0.y, W1.z) + Dither * DitherWeights[3] * GlobalDitherNoiseWeight)
		   + SampleLPVColor(vec3(W0.x, W0.y, W1.z) - Dither * DitherWeights[3] * GlobalDitherNoiseWeight);
	return max((Interpolated / 8.0), 0.00000001f);
}


float InterpLPVDensity(vec3 UV) 
{
    vec3 LPVResolution = vec3(384.0f, 128.0f, 384.0f);
    vec3 FractTexel = fract(UV * LPVResolution);
    vec3 LinearOffset = (FractTexel * (FractTexel - 1.0f) + 0.5f) / LPVResolution;
    vec3 W0 = UV - LinearOffset;
    vec3 W1 = UV + LinearOffset;
    float Density = texture(u_LPVLightLevel, vec3(W0.x, W0.y, W0.z)).x
    	          + texture(u_LPVLightLevel, vec3(W1.x, W0.y, W0.z)).x
    	          + texture(u_LPVLightLevel, vec3(W1.x, W1.y, W0.z)).x
    	          + texture(u_LPVLightLevel, vec3(W0.x, W1.y, W0.z)).x
    	          + texture(u_LPVLightLevel, vec3(W0.x, W1.y, W1.z)).x
    	          + texture(u_LPVLightLevel, vec3(W1.x, W1.y, W1.z)).x
    	          + texture(u_LPVLightLevel, vec3(W1.x, W0.y, W1.z)).x
		          + texture(u_LPVLightLevel, vec3(W0.x, W0.y, W1.z)).x;
	return max(Density / 8.0, 0.00000001f);
}



vec3 GetSmoothLPVData(vec3 UV) {    
    UV *= 1.0f/vec3(384.0f,128.0f,384.0f);
    //return vec3(InterpLPVDensity(UV)*50.0f)*pow(InterpolateLPVColorData(UV),vec3(1.0f/1.8f))*2.0f;
    return vec3(InterpLPVDensity(UV)*50.0f)*pow(InterpolateLPVColorDithered(UV),vec3(1.0f/1.8f))*2.0f;
}

vec3 GetSmoothLPVDensity(vec3 UV) {    
    UV *= 1.0f/vec3(384.0f,128.0f,384.0f);
    return vec3(InterpLPVDensity(UV)*54.0f);
}


float SRGBToLinear(float x){
    return x > 0.04045 ? pow(x * (1 / 1.055) + 0.0521327, 2.4) : x / 12.92;
}

vec3 SRGBToLinearVec3(vec3 x){
    return vec3(SRGBToLinear(x.x),
                SRGBToLinear(x.y),
                SRGBToLinear(x.z));
}

vec3 TemperatureToRGB(float temperatureInKelvins)
{
	vec3 retColor;
	
    temperatureInKelvins = clamp(temperatureInKelvins, 1000, 50000) / 100;
    
    if (temperatureInKelvins <= 66){
        retColor.r = 1;
        retColor.g = clamp01(0.39008157876901960784 * log(temperatureInKelvins) - 0.63184144378862745098);
    } else {
    	float t = temperatureInKelvins - 60;
        retColor.r = clamp01(1.29293618606274509804 * pow(t, -0.1332047592));
        retColor.g = clamp01(1.12989086089529411765 * pow(t, -0.0755148492));
    }
    
    if (temperatureInKelvins >= 66)
        retColor.b = 1;
    else if(temperatureInKelvins <= 19)
        retColor.b = 0;
    else
        retColor.b = clamp01(0.54320678911019607843 * log(temperatureInKelvins - 10) - 1.19625408914);

    return SRGBToLinearVec3(retColor);
}     



vec2 ProjectDirection(vec3 Direction, vec2 TextureSize) 
{
    //vec2 TextureSize = vec2(textureSize(u_DebugTexture,0).xy);
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

// end of light combine / color pass