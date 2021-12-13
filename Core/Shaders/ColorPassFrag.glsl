// Light combine / Color pass.

#version 450 core

#define CLOUD_HEIGHT 70
#define PI 3.14159265359
#define THRESH 1.41414

#define clamp01(x) (clamp(x, 0.0f, 1.0f))

// Outputs 
layout (location = 0) out vec3 o_Color;
layout (location = 1) out vec4 o_PBR;

// Vertex shader inputs 
vec2 g_TexCoords;
in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

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


uniform sampler2D u_SSSShadowMap;
uniform bool u_SSSSS;
uniform float u_SubsurfaceScatterStrength;



uniform sampler2D u_DiffuseSHy;
uniform sampler2D u_DiffuseCoCg;

uniform float u_TextureDesatAmount;

uniform sampler2D u_ReflectionSHData;
uniform sampler2D u_ReflectionCoCgData;

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
uniform int u_GrassBlockProps[10];

uniform vec2 u_Halton;


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
vec4 textureBicubic(sampler2D sampler, vec2 texCoords);
vec3 CalculateDirectionalLight(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 radiance_s, vec3 albedo, vec3 normal, vec3 pbr, float shadow);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
vec3 GetSmoothLPVData(vec3 UV);
vec3 GetSmoothLPVDensity(vec3 UV);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv, out vec4 UVDerivative);
bool CompareFloatNormal(float x, float y);
vec3 GetNormalFromID(float n);
vec3 SampleNormal(sampler2D samp, vec2 txc);
vec3 GetRayDirectionAt(vec2 screenspace);
vec4 SamplePositionAt(sampler2D pos_tex, vec2 txc);
vec4 GetTextureIDs(int BlockID);
int GetBlockID(vec2 txc);
bool InScreenSpace(vec2 x);
vec2 hash2();
vec3 SampleCone(vec2 Xi, float CosThetaMax); 
vec3 XYZToRGB(in vec3 xyz);
vec3 RGBToXYZ(in vec3 rgb);
float MiePhaseFunction(float x, float g);
vec3 fresnelroughness(vec3 Eye, vec3 norm, vec3 F0, float roughness);
float FastDistance(vec3 p1, vec3 p2);
vec4 BetterTexture(sampler2D samp, vec2 uv);
vec3 BasicSaturation(vec3 Color, float Adjustment);
float RayCapsuleIntersection(in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float r);
vec3 saturate(vec3 x);
vec2 SmoothStepUV(vec2 uv, vec2 res, float width);
float GetLuminance(vec3 color);
vec4 ClampedTexture(sampler2D tex, vec2 txc);
vec4 texture_catmullrom(sampler2D tex, vec2 uv);





// Constants / Globals 


const vec3 ATMOSPHERE_SUN_COLOR = vec3(1.0f);
const vec3 ATMOSPHERE_MOON_COLOR =  vec3(0.1f, 0.1f, 1.0f);
const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);
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

    vec3 ray_dir = normalize(in_ray_dir);
    
    if (true_transmittance > 0.325f) {

        if(dot(ray_dir, sun_dir) > 0.999825f)
        {
            atmosphere_color = ATMOSPHERE_SUN_COLOR; 
            o_PBR.w = float(1.0f);
            return true;
        }
        
        if(dot(ray_dir, moon_dir) > 0.99986f)
        {
            atmosphere_color = ATMOSPHERE_MOON_COLOR;
            o_PBR.w = float(1.2f);
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
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;

    // Calculate colors ->
    vec3 BaseSun = SAMPLED_SUN_COLOR * vec3(1.0f, 0.9f, 0.9f)* 1.0f;
    float DuskVisibility = clamp(pow(abs(u_SunDirection.y - 1.0), 2.0f), 0.0f, 1.0f);
    BaseSun = mix(vec3(2.2f), BaseSun, DuskVisibility);
    vec3 ScatterColor = mix(BaseSun, BasicSaturation(SAMPLED_MOON_COLOR_RAW, 1.64f)*1.6*(0.9f/1.0f), SunVisibility); 

    // Normalize ->
    ScatterColor *= 1.0f / PI;
    ScatterColor *= 1.2f;
    ScatterColor = clamp(ScatterColor, 0.0f, 2.0f);


    vec3 NormalizedDir = normalize(v_RayDirection);
    float Dither = bayer128(gl_FragCoord.xy);

    vec2 SampleCoord = g_TexCoords + (Dither / 768.0f);


	vec4 SampledCloudData;

    if (u_CloudCatmullRomUpsampling) {
        SampledCloudData = texture_catmullrom(u_CloudData, SampleCoord).rgba;
	}
    
    else {
        SampledCloudData = BetterTexture(u_CloudData, SampleCoord).rgba;
    }

    //vec4 Detail = UpscaleCloudDetail(SampledCloudData);
    //SampledCloudData.xyz = mix(SampledCloudData.xyz, Detail.xyz, 1.0f - SampledCloudData.w);

    vec3 Sky = vec3(0.0f);

    bool v = GetAtmosphere(Sky, NormalizedDir, SampledCloudData.w*20.5f, SampledCloudData.w);

    if (!u_CloudsEnabled) {

        return Sky;
    }

    float transmittance = max(SampledCloudData.w, 0.07525f);
    return (Sky * transmittance) + (((SampledCloudData.xyz * ScatterColor) + (Dither / 1024.0f)));

}

vec3 FetchSky() { return GetAtmosphereAndClouds(); }

vec2 ReprojectShadow(in vec3 world_pos)
{
	vec3 WorldPos = world_pos;

	vec4 ProjectedPosition = u_ShadowProjection * u_ShadowView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}


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

    vec3 BiasedWorldPos = world_pos - (flat_normal * 0.3f);
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

//
float GetDisplacementAt(in vec2 txc, in float pbridx) 
{
    return texture(u_BlockPBRTextures, vec3(vec2(txc.x, txc.y), pbridx)).b * 0.35f * u_POMHeight;
}

// Parallax occlusion mapping 
// Exponential ray step 
// Tried dithering ray step, which resulted in not-good results (introduces noise which is hard to tackle as you need this to be temporally coherent)
// - to avoid artifacts


vec2 ParallaxOcclusionMapping(vec2 TextureCoords, vec3 ViewDirection, in float pbridx) // View direction should be in tangent space!
{ 
    if (u_DitherPOM) {
        float Bayer = bayer64(gl_FragCoord.xy);
        ViewDirection *= mix(0.9f, 1.0f, Bayer);
    }

    float NumLayers = u_HighQualityPOM ? 96 : 64; 
    float LayerDepth = 1.0 / (NumLayers * 0.65);
    float CurrentLayerDepth = 0.0;
    vec2 P = ViewDirection.xy * 1.0f; 
    vec2 DeltaTexCoords = P / NumLayers;
    vec2 InitialDeltaCoords = DeltaTexCoords;

    vec2  CurrentTexCoords = TextureCoords;
    float CurrentDepthMapValue = GetDisplacementAt(CurrentTexCoords, pbridx);

    for (int i = 0 ; i < NumLayers ; i++)
    {
        if(CurrentLayerDepth < CurrentDepthMapValue)
        {
            //CurrentTexCoords -= max(DeltaTexCoords, InitialDeltaCoords * 0.025f);
            CurrentTexCoords -= DeltaTexCoords;
            CurrentDepthMapValue = GetDisplacementAt(CurrentTexCoords, pbridx);  
            CurrentLayerDepth += LayerDepth;
            DeltaTexCoords *=  0.95f;
        }
    }

    vec2 PrevTexCoords = CurrentTexCoords + DeltaTexCoords;
    float AfterDepth  = CurrentDepthMapValue - CurrentLayerDepth;
    float BeforeDepth = GetDisplacementAt(PrevTexCoords, pbridx) - CurrentLayerDepth + LayerDepth;
    float Weight = AfterDepth / (AfterDepth - BeforeDepth);
    vec2 FinalTexCoords = PrevTexCoords * Weight + CurrentTexCoords * (1.0 - Weight);
    return FinalTexCoords;
}   




bool IsAtEdge(in vec2 txc)
{
    vec2 TexelSize = 1.0f / textureSize(u_InitialTracePositionTexture, 0);

    const vec2 Kernel[8] = vec2[8]
	(
		vec2(-1.0, -1.0),
		vec2( 0.0, -1.0),
		vec2( 1.0, -1.0),
		vec2(-1.0,  0.0),
		vec2( 1.0,  0.0),
		vec2(-1.0,  1.0),
		vec2( 0.0,  1.0),
		vec2( 1.0,  1.0)
	);

    for (int i = 0 ; i < 8 ; i++)
    {
        if (SamplePositionAt(u_InitialTracePositionTexture, txc + Kernel[i] * TexelSize).w <= 0.0f)
        {
            return true;
        }
    }

    return false;
}

float SHToY(vec4 shY)
{
    return max(0, 3.544905f * shY.w); // get luminance (Y) from the first spherical harmonic
}

// Spatially upscales indirect data ->
void SpatialUpscaleData(vec3 BaseNormal, float BaseLinearDepth, out vec4 SH, out vec2 CoCg, out vec4 SpecularSHy, out vec2 SpecularCoCg, out float ShadowSample, bool fuckingsmooth)
{
    const bool BE_FUCKING_USELESS = false;

    if (BE_FUCKING_USELESS) {
        vec2 SampleCoord = g_TexCoords;
        SH = texture(u_DiffuseSHy, SampleCoord).xyzw;
		CoCg += texture(u_DiffuseCoCg, SampleCoord).xy;
        SpecularSHy += texture(u_ReflectionSHData, SampleCoord).xyzw;
        SpecularCoCg += texture(u_ReflectionCoCgData, SampleCoord).xy;
        ShadowSample += texture(u_ShadowTexture, SampleCoord).x;
        return;
    }

	vec4 TotalSH = vec4(0.0f);
	vec2 TotalCoCg = vec2(0.0f);
    vec4 TotalSpecSH = vec4(0.0f);
    vec2 TotalSpecCoCg = vec2(0.0f);
	float TotalWeight = 0.0f;
	//float TotalShadowWeight = 0.0f;

    const vec2 Offsets[5] = vec2[5](
        vec2(0.0f, 1.0f),
        vec2(1.0f, 0.0f),
        vec2(0.0f),
        vec2(-1.0f, 0.0f),
        vec2(0.0, -1.0f)
    );

    //vec4 SpecularSamplesSHy[5];
    //vec2 SpecularSamplesCoCg[5];

	vec2 TexelSize = 1.0f / textureSize(u_DiffuseSHy, 0);
	const float Weights[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);
	//float BaseShadowSample = texture(u_ShadowTexture, g_TexCoords).x;
	//float BaseL = GetLuminance(vec3(BaseShadowSample));

	for (int s = 0 ; s < 5 ; s++) {
		
        vec2 SampleCoord = g_TexCoords + (vec2(Offsets[s])) * TexelSize;
		float LinearDepthAt = (texture(u_InitialTracePositionTexture, SampleCoord).x);
        
        float DepthWeight = 1.0f / (abs(LinearDepthAt - BaseLinearDepth) + 0.01f);
		DepthWeight = pow(DepthWeight, 16.0f);
        DepthWeight = max(DepthWeight, 0.0f);

        vec3 NormalAt = SampleNormal(u_NormalTexture, SampleCoord.xy).xyz;
		float NormalWeight = pow(max(abs(dot(NormalAt, BaseNormal)),0.000001f), 64.0f);
        NormalWeight = max(NormalWeight, 0.0f);

		float Weight = clamp(NormalWeight * DepthWeight, 0.0f, 1.0f);

		Weight = max(Weight, 0.01f);
        float SpecularWeight = Weight;

		TotalSH += texture(u_DiffuseSHy, SampleCoord).xyzw * Weight;
		TotalCoCg += texture(u_DiffuseCoCg, SampleCoord).xy * Weight;
        //SpecularSamplesSHy[s] = texture(u_ReflectionSHData, SampleCoord).xyzw;
        //SpecularSamplesCoCg[s] = texture(u_ReflectionCoCgData, SampleCoord).xy;

        TotalSpecSH += texture(u_ReflectionSHData, SampleCoord).xyzw * Weight;
        TotalSpecCoCg += texture(u_ReflectionCoCgData, SampleCoord).xy * Weight;
		
		//float ShadowSampleAt = texture(u_ShadowTexture, SampleCoord).x;
		//float ShadowWeight = Weight * clamp(pow(abs(BaseL-GetLuminance(vec3(ShadowSampleAt)))/2.0f,64.0f), 0.025f, 1.0f);
        //ShadowSample += ShadowSampleAt * ShadowWeight;

		TotalWeight += Weight;
		//TotalShadowWeight += ShadowWeight;
	}

    TotalWeight = max(TotalWeight, 0.01f);

    SH = TotalSH / TotalWeight;
    CoCg = TotalCoCg / TotalWeight;
    //ShadowSample = ShadowSample / TotalShadowWeight;
	ShadowSample = BetterTexture(u_ShadowTexture,g_TexCoords).x;
    TotalSpecSH = TotalSpecSH / TotalWeight;
    TotalSpecCoCg = TotalSpecCoCg / TotalWeight;
    SpecularSHy = TotalSpecSH;
    SpecularCoCg = TotalSpecCoCg;
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

// Approximates dfg term ->
vec3 DFGPolynomialApproximate(vec3 F0, float NdotV, float roughness) // from unreal 
{
    const vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
    const vec4 c1 = vec4(1.0, 0.0425, 1.04, -0.04);
    vec4 r = roughness * c0 + c1;
    float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
    vec2 AB = vec2(-1.04, 1.04) * a004 + r.zw;
    return F0 * AB.x + AB.y;
}


// cook torrance brdf : 

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
	vec3 F = fresnelroughness(Lo, normal.xyz, vec3(F0), pbr.r); // use fresnel roughness, gives a slightly better approximation
	float D = ndfGGX(cosLh, pbr.r);
	float G = gaSchlickGGX(cosLi, cosLo, pbr.r);

	vec3 kd = mix(vec3(1.0) - F, vec3(0.0), pbr.g);
	vec3 diffuseBRDF = kd * albedo;

	vec3 specularBRDF = (F * D * G) / max(Epsilon, 4.0 * cosLi * cosLo); specularBRDF = clamp(specularBRDF, 0.0f, 2.0f);
	vec3 Result = (diffuseBRDF * Lradiance * cosLi) + (specularBRDF * radiance_s * cosLi);
    return clamp(Result, 0.0f, 2.5) * clamp((1.0f - Shadow), 0.0f, 1.0f);
}


vec3 TemperatureToRGB(float temperatureInKelvins);

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
    float TimeOfDayTransition = sv * sv;
    vec3 MoonColor = MoonTransmittance;

    if (u_NebulaCelestialColor)
    {
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
    MoonColor = mix(MoonColor, vec3(L), 0.05f); 
    return MoonColor * 0.5f;
}

// Custon very non physically based subsurface scattering 
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


vec2 BetterFilteringUVTweak(vec2 res, vec2 uv) 
{
    uv = uv*res + 0.5;
    vec2 fl = floor(uv);
    vec2 fr = fract(uv);
    vec2 aa = fwidth(uv)*0.75;
    fr = smoothstep( vec2(0.5)-aa, vec2(0.5)+aa, fr);
    uv = (fl+fr-0.5) / res;
    return uv;
}

void main()
{
    g_TexCoords = v_TexCoords;
    //g_TexCoords += u_Halton/u_Dimensions;


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

    if (WorldPosition.w > 0.0f)
    {
       

        float Bias = 0.0035f;
        bool InBiasedSS =  (g_TexCoords.x > 0.0 + Bias && g_TexCoords.x < 1.0 - Bias 
		 && g_TexCoords.y > 0.0 + Bias && g_TexCoords.y < 1.0f - Bias);

        //if (!IsAtEdge(g_TexCoords))
        if (true)
        {
            vec2 UV;
            vec3 Tangent, Bitangent;
            vec4 UVDerivative;

            CalculateVectors(WorldPosition.xyz, SampledNormals, Tangent, Bitangent, UV, UVDerivative); 

	        mat3 tbn = mat3((Tangent), (Bitangent), (SampledNormals));
            int BaseBlockID = GetBlockID(g_TexCoords);
            bool DoSSS = BlockSSSSSData[BaseBlockID] > 0 && u_SSSSS;
            vec4 data = GetTextureIDs(BaseBlockID);


            // Handle grass block! 
            // dont kill me pls :( 
            if (BaseBlockID == u_GrassBlockProps[0])
	        {
	            if (SampledNormals == NORMAL_LEFT || SampledNormals == NORMAL_RIGHT || SampledNormals == NORMAL_FRONT || SampledNormals == NORMAL_BACK)
	        	{
	        		data.x = u_GrassBlockProps[4];
	        		data.y = u_GrassBlockProps[5];
	        		data.z = u_GrassBlockProps[6];
	        	}

	        	else if (SampledNormals == NORMAL_TOP)
	        	{
	        		data.x = u_GrassBlockProps[1];
	        		data.y = u_GrassBlockProps[2];
	        		data.z = u_GrassBlockProps[3];
	        	}

	        	else if (SampledNormals == NORMAL_BOTTOM)
	        	{
	        		data.x = u_GrassBlockProps[7];
	        		data.y = u_GrassBlockProps[8];
	        		data.z = u_GrassBlockProps[9];
	        	}
	        }



            if (u_POM && data.r != u_GrassblockAlbedoID)
            {
                vec2 InitialUV = UV;

                if (SampledNormals == vec3(-1.0f, 0.0f, 0.0f)) {  UV.x = 1.0f - UV.x; UV.y = 1.0f - UV.y; }
                if (SampledNormals == vec3(1.0f, 0.0f, 0.0f)) {  UV.y = 1.0f - UV.y; }
                if (SampledNormals == vec3(0.0f, -1.0f, 0.0f)) {  UV.y = 1.0f - UV.y; }
                //if (SampledNormals == vec3(0.0f, 0.0f, 1.0f)) {  flip = true;}
                //if (SampledNormals == vec3(0.0f, 0.0f, -1.0f)) {  flip = true; }

                vec3 TangentViewPosition = tbn * u_ViewerPosition;
                vec3 TangentFragPosition = tbn * WorldPosition.xyz; 
                vec3 TangentViewDirection = normalize((TangentFragPosition - TangentViewPosition));
                
                UV = ParallaxOcclusionMapping(UV, TangentViewDirection, data.z);
                //UV.y = flip && data.r == u_GrassblockAlbedoID ? 1.0f - UV.y : UV.y;
            
            } else { UV.y = 1.0f - UV.y; }

            //vec3 AlbedoColor = texture(u_BlockAlbedoTextures, vec3(UV.xy, data.x)).rgb;
            //vec4 PBRMap = texture(u_BlockPBRTextures, vec3(UV, data.z)).rgba;
            //vec3 NormalMapped = tbn * (texture(u_BlockNormalTextures, vec3(UV, data.y)).rgb * 2.0f - 1.0f);

            vec2 SmoothstepUV = BetterFilteringUVTweak(vec2(512.0f), UV);

            // Fix texture seam with approximated derivative ->
            
            vec4 PBRMap = textureGrad(u_BlockPBRTextures, vec3(SmoothstepUV, data.z), UVDerivative.xy, UVDerivative.zw).rgba;
            vec3 AlbedoColor = textureGrad(u_BlockAlbedoTextures, vec3(SmoothstepUV, data.x), UVDerivative.xy, UVDerivative.zw).rgb;
            vec3 NormalMapped = tbn * (textureGrad(u_BlockNormalTextures, vec3(SmoothstepUV, data.y), UVDerivative.xy, UVDerivative.zw).rgb * 2.0f - 1.0f);
            
            AlbedoColor = BasicSaturation(AlbedoColor, 1.0f - u_TextureDesatAmount);

            if (PBRMap.y > 0.02f){
                AlbedoColor = BasicSaturation(AlbedoColor, 0.8f);
            }

            vec3 NonAmplifiedNormal = NormalMapped;

            if (u_AmplifyNormalMap) {
                NormalMapped.x *= 1.64f;
                NormalMapped.z *= 1.85f;
                NormalMapped += 1e-4f;
                NormalMapped = normalize(NormalMapped);
            }
            
            float Emissivity = data.w > -0.5f ? textureLod(u_BlockEmissiveTextures, vec3(UV, data.w), 0.0f).r : 0.0f;

            //vec4 Diffuse = BilateralUpsample(u_DiffuseTexture, g_TexCoords, SampledNormals.xyz, WorldPosition.z);
            //vec4 Diffuse = PositionOnlyBilateralUpsample(u_DiffuseTexture, g_TexCoords, WorldPosition.xyz);
            //vec3 Diffuse = DepthOnlyBilateralUpsample(u_DiffuseTexture, g_TexCoords, WorldPosition.z).xyz;
            //vec4 SampledIndirectDiffuse = BilateralUpsample2(u_DiffuseTexture, g_TexCoords, WorldPosition.xyz, SampledNormals.xyz).xyzw;
            
            vec4 SHy;
            vec2 ShCoCg;
            vec4 SpecularSH;
            vec2 SpecularCoCg;
			float UpscaledShadow=0.0f;

		
            SpatialUpscaleData(SampledNormals.xyz, WorldPosition.w, SHy, ShCoCg, SpecularSH, SpecularCoCg,UpscaledShadow, PBRMap.x <= 0.1);
			UpscaledShadow = clamp(UpscaledShadow,0.0f,1.0f);

			float RayTracedShadow = ComputeShadow(WorldPosition.xyz, SampledNormals.xyz,UpscaledShadow);
			RayTracedShadow = clamp(RayTracedShadow,0.,1.);


            vec3 IndirectN = NormalMapped.xyz;
            vec3 SampledIndirectDiffuse = SHToIrridiance(SHy, ShCoCg, IndirectN);

            // VXAO from indirect diffuse trace : 
            // 

            bool do_vxao = u_DoVXAO && u_SVGFEnabled;
			
			const int VXGI_CUTOFF = 69 + 7;
			
            if ((do_vxao && !u_RTAO && (distance(WorldPosition.xyz, u_ViewerPosition) < VXGI_CUTOFF)) || (!u_VXAOCutoff)) // -> Causes artifacts if the AO is applied too far away
            {
                float bias = 0.00125f;
                if (g_TexCoords.x > bias && g_TexCoords.x < 1.0f - bias &&
                    g_TexCoords.y > bias && g_TexCoords.y < 1.0f - bias)
                {
                    float ao = BetterTexture(u_VXAO, g_TexCoords).x;
					
					// Multiply indirect diffuse by ao 
                    SampledIndirectDiffuse.xyz *= vec3(clamp(pow(ao, 1.85f), 0.1f, 1.0f));
                }
            }

            vec3 LightAmbience = (vec3(120.0f, 172.0f, 255.0f) / 255.0f) * 1.01f;
            vec3 Ambient = (AlbedoColor * LightAmbience) * 0.09f;
            float SampledAO = pow(PBRMap.w, 1.25f); // Ambient occlusion map

            // Interpolate and find sun colors : 
            //float DuskVisibility = clamp(pow(distance(u_SunDirection.y, 1.0), 1.5f), 0.0f, 1.0f);

             //mix(SUN_COLOR, DUSK_COLOR * 0.5f, DuskVisibility);

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
            

            if (InBiasedSS) {
                vec3 I = normalize(WorldPosition.xyz - u_ViewerPosition);
                vec3 R = reflect(I, NormalMapped).xyz;
                SpecularIndirect += SHToIrridiance(SpecularSH, SpecularCoCg, normalize(R)); 
            }

           
			SpecularIndirect = max(vec3(0.0001f), SpecularIndirect);
            vec3 Lo = normalize(u_ViewerPosition - WorldPosition.xyz); // Outgoing direction 
            vec3 F0 = mix(vec3(0.04), AlbedoColor, PBRMap.g); // Fresnel at 0 degrees.

            vec3 SpecularFactor = fresnelroughness(Lo, NormalMapped.xyz, vec3(F0), Roughness); 
            
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

           //o_Color = vec3(SubsurfaceScatter);
            //o_Color = vec3(NormalMapped);
           // return;

            if (!dfg) {
                o_Color = ((DirectLighting + SubsurfaceScatter) + ((1.0f - SpecularFactor) * DiffuseIndirect) + 
                      (SpecularFactor * SpecularIndirect));
            }

            else {

                // direct + ind_diff + ind_spec * dfg (fr) 
                // dfg * fresnel to correct fresnel for energy conservation
                // 1.0 + f0 * (1.0 / dfg.y - 1.0) usually

                float NDotV = clamp(dot(NonAmplifiedNormal, Lo), 0.0f, 1.0f);
                vec3 Tint = mix(vec3(1.0f), AlbedoColor, float(PBRMap.g > 0.0125f));
                vec3 DFG = DFGPolynomialApproximate(F0, Roughness, NDotV); // dfg
                vec3 IndirectSpecularBRDF = DFG * clamp(Tint, 0.0f, 1.0f);

                o_Color = (DirectLighting + SubsurfaceScatter) + DiffuseIndirect * (1.0f - IndirectSpecularBRDF);
                o_Color += SpecularIndirect * IndirectSpecularBRDF;
            }

            
			

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
          
            // Output utility : 
            //o_Normal = vec3(NormalMapped.x, NormalMapped.y, NormalMapped.z);
            o_PBR.xyz = PBRMap.xyz;
            o_PBR.w = Emissivity;


            // Fix bloom light leak around the edges : 
            const bool BloomLightLeakFix = true;
            if (BloomLightLeakFix) {
                float lbiasx = 0.032501f;
                float lbiasy = 0.032501f;
                o_PBR.w *= float(UV.x > lbiasx && UV.x < 1.0f - lbiasx &&
                                 UV.y > lbiasy && UV.y < 1.0f - lbiasy);
            }

            //return;
        }

        else 
        {   
            vec3 CloudAndSky = GetAtmosphereAndClouds();
			o_Color = (CloudAndSky);
			//o_Normal = vec3(-1.0f);
			o_PBR.xyz = vec3(-1.0f);
        }
    }

    else 
    {   
        vec3 CloudAndSky = GetAtmosphereAndClouds();
        o_Color = (CloudAndSky);
        //o_Normal = vec3(-1.0f);
        o_PBR.xyz = vec3(-1.0f);
    }

    if (!true) {
        //vec4 DebugData = texture(u_DebugTexture, ProjectDirection(v_RayDirection)).xyzw;
        vec4 DebugData = texture(u_DebugTexture, g_TexCoords).xyzw;
        //o_Color = max(vec3(DebugData.x / 20.0f), 0.0f);
        o_Color = max(vec3(DebugData.x)/20.0,0.00001f);
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
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

vec4 GetTextureIDs(int BlockID) 
{
	return vec4(float(BlockAlbedoData[BlockID]),
				float(BlockNormalData[BlockID]),
				float(BlockPBRData[BlockID]),
				float(BlockEmissiveData[BlockID]));
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

vec3 fresnelroughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
{
    float cosTheta = max(dot(norm, Eye), 0.0);
    const float magic = 2.4f;
    return F0 + (max(vec3(pow(1.0f - roughness, magic)), F0) - F0) * pow(clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
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

vec2 SmoothStepUV(vec2 uv, vec2 res, float width)
{
    uv = uv * res;
    vec2 uv_floor = floor(uv + 0.5);
    vec2 uv_fract = fract(uv + 0.5);
    vec2 uv_aa = fwidth(uv) * width * 0.5;
    uv_fract = smoothstep(vec2(0.5) - uv_aa, vec2(0.5) + uv_aa, uv_fract);
    return (uv_floor + uv_fract - 0.5) / res;
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

vec3 InterpolateLPVColorDataBicubic(vec3 position) 
{
    position *= vec3(384.0f, 128.0f, 384.0f);
    position = position - 0.5;

	vec3 f = fract(position);
	position -= f;

	vec3 ff = f * f;
	vec3 w0_1;
	vec3 w0_2;
	vec3 w1_1;
	vec3 w1_2;

	w0_1 = 1 - f; w0_1 *= w0_1 * w0_1;
	w1_2 = ff * f;
	w1_1 = 3 * w1_2 + 4 - 6 * ff;
	w0_2 = 6 - w1_1 - w1_2 - w0_1;

	vec3 s1 = w0_1 + w1_1;
	vec3 s2 = w0_2 + w1_2;
	vec3 cLo = position.xyz + vec3(-0.5f) + w1_1 / s1;
	vec3 cHi = position.xyz + vec3(1.5f) + w1_2 / s2;
	vec3 m = s1 / (s1 + s2);
	m = 1.0 - m;

	vec3 lightLoYLoZ = mix(SampleLPVColorTexel(ivec3(cLo.x, cLo.y, cLo.z)), SampleLPVColorTexel(ivec3(cHi.x, cLo.y, cLo.z)), m.x);
	vec3 lightLoYHiZ = mix(SampleLPVColorTexel(ivec3(cLo.x, cLo.y, cHi.z)), SampleLPVColorTexel(ivec3(cHi.x, cLo.y, cHi.z)), m.x);
	vec3 lightHiYLoZ = mix(SampleLPVColorTexel(ivec3(cLo.x, cHi.y, cLo.z)), SampleLPVColorTexel(ivec3(cHi.x, cHi.y, cLo.z)), m.x);
	vec3 lightHiYHiZ = mix(SampleLPVColorTexel(ivec3(cLo.x, cHi.y, cHi.z)), SampleLPVColorTexel(ivec3(cHi.x, cHi.y, cHi.z)), m.x);
	vec3 lightLoZ = mix(lightLoYLoZ, lightHiYLoZ, m.y);
	vec3 lightHiZ = mix(lightLoYHiZ, lightHiYHiZ, m.y);
	return mix(lightLoZ, lightHiZ, m.z);
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
    const float DitherWeights[4] = float[4](1.0f, 0.5f, 0.35f, 0.25f);
    const float GlobalDitherNoiseWeight = 3.0f;

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
    return vec3(InterpLPVDensity(UV)*50.0f)*pow(InterpolateLPVColorData(UV),vec3(1.0f/1.8f))*2.0f;
    //return vec3(InterpLPVDensity(UV)*50.0f)*pow(InterpolateLPVColorDithered(UV),vec3(1.0f/1.8f))*2.0f;
}

vec3 GetSmoothLPVDensity(vec3 UV) {    
    UV *= 1.0f/vec3(384.0f,128.0f,384.0f);
    return vec3(InterpLPVDensity(UV)*54.0f);
}

// Bicubic

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

// catmull rom interpolation.
#define sqr(x) (x*x)

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

bool CompareVec3(vec3 v1, vec3 v2) {
	float e = 0.0125f;
	return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
}

// Use a derivative to fix the uv seam ->
vec4 GetUVDerivative(in vec2 P) 
{
    vec2 UV = fract(P);
    vec2 dx = dFdx(UV);
    vec2 dy = dFdy(UV);
    vec2 OffsettedUV_1  = fract(UV + 0.25f);
	vec2 OffsettedUV_2 = fract(UV + 0.5f);
    vec2 dx2 = dFdx(OffsettedUV_1);
	vec2 dy2 = dFdy(OffsettedUV_1);

	if(dot(dx, dx) + dot(dy, dy) > dot(dx2, dx2) + dot(dy2, dy2)) 
    {
		dx = dx2;
		dy = dy2;
	}

    return vec4(dx, dy);
}

void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );

	const vec3 Tangents[6] = vec3[]( vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(0.0f, 0.0f, -1.0f), vec3(0.0f, 0.0f, -1.0f)
				   );

	const vec3 BiTangents[6] = vec3[]( vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f),
				     vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 1.0f),
					 vec3(0.0f, -1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f)
	);

    uv = vec2(1.0f);

	if (CompareVec3(normal, Normals[0]))
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[0];
		bitangent = BiTangents[0];
    }

    else if (CompareVec3(normal, Normals[1]))
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[1];
		bitangent = BiTangents[1];
    }

    else if (CompareVec3(normal, Normals[2]))
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[2];
		bitangent = BiTangents[2];
    }

    else if (CompareVec3(normal, Normals[3]))
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[3];
		bitangent = BiTangents[3];
    }
	
    else if (CompareVec3(normal, Normals[4]))
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[4];
		bitangent = BiTangents[4];
    }
    

    else if (CompareVec3(normal, Normals[5]))
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[5];
		bitangent = BiTangents[5];
    }
}

void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv, out vec4 UVDerivative)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );

	const vec3 Tangents[6] = vec3[]( vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(0.0f, 0.0f, -1.0f), vec3(0.0f, 0.0f, -1.0f)
				   );

	const vec3 BiTangents[6] = vec3[]( vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f),
				     vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 1.0f),
					 vec3(0.0f, -1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f)
	);

    uv = vec2(1.0f);

	if (CompareVec3(normal, Normals[0]))
    {
        uv = vec2(fract(world_pos.xy));
        UVDerivative = GetUVDerivative(world_pos.xy);
		tangent = Tangents[0];
		bitangent = BiTangents[0];
    }

    else if (CompareVec3(normal, Normals[1]))
    {
        uv = vec2(fract(world_pos.xy));
        UVDerivative = vec4(GetUVDerivative(world_pos.xy));
		tangent = Tangents[1];
		bitangent = BiTangents[1];
    }

    else if (CompareVec3(normal, Normals[2]))
    {
        uv = vec2(fract(world_pos.xz));
        UVDerivative = vec4(GetUVDerivative(world_pos.xz));
		tangent = Tangents[2];
		bitangent = BiTangents[2];
    }

    else if (CompareVec3(normal, Normals[3]))
    {
        uv = vec2(fract(world_pos.xz));
        UVDerivative = vec4(GetUVDerivative(world_pos.xz));
		tangent = Tangents[3];
		bitangent = BiTangents[3];
    }
	
    else if (CompareVec3(normal, Normals[4]))
    {
        uv = vec2(fract(world_pos.zy));
        UVDerivative = vec4(GetUVDerivative(world_pos.zy));
		tangent = Tangents[4];
		bitangent = BiTangents[4];
    }
    

    else if (CompareVec3(normal, Normals[5]))
    {
        uv = vec2(fract(world_pos.zy));
        UVDerivative = vec4(GetUVDerivative(world_pos.zy));
		tangent = Tangents[5];
		bitangent = BiTangents[5];
    }
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