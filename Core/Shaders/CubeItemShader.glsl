#version 430 core 
#define INF 100000.0f
#define PI 3.14159265359f
#define TAU (2.0f * PI)
#define clamp01(x) (clamp(x,0.,1.))

// Shape :
// The cube is ray casted using a simple in-shader camera and rotation matrix

// Lighting :
// Lighting is done completely stylistically 
// It uses a simple lambertian BRDF for the direct lighting 
// Indirect lighting is done using a cut down version of image based lighting 
// Indirect diffuse uses a pre convoluted irradiance map 
// Specular is calculated in realtime by sampling the ggx vndf and using a bunch of samples.

// Antialiasing :
// Uses a custom anti aliasing method 
// -> Use screen space derivatives to find out luminance jumps
// And only apply anti aliasing where it jumps. 
// The anti aliasing algorithm samples multiple subpixels linearly from -0.5 -> +0.5

// I spent around 2 days working on this shit. What am I doing with my fucking life.


//#define SMART_ANTI_ALIAS_DEBUG


layout(rgba16f, binding = 1) uniform image2D o_OutputImage;

// Screen space uv ->
in vec2 v_TexCoords;

uniform float u_Time;
uniform vec2 u_Dimensions;
uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;

// Sky/HDRI ->
uniform samplerCube u_Sky;
uniform samplerCube u_RandomHDRI;
uniform samplerCube u_RandomHDRIDiffuse;

// Block textures ->
uniform sampler2DArray u_AlbedoTextures;
uniform sampler2DArray u_NormalTextures;
uniform sampler2DArray u_EmissiveTextures;
uniform sampler2DArray u_PBRTextures;

uniform int u_HeldBlockID;

// Params
uniform int u_AntialiasLevel;
uniform int u_SpecularSampleBias;
uniform bool u_SimpleLighting;

uniform float u_SunVisibility;

// Multitexturing *FUCKING HACK* for only the grass block
// **Will be removed in the future**
uniform int u_GrassBlockProps[10];


// IDs ->
layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
	int BlockSSSSSData[128];
};

// Constants 
const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

// Calculates tangent, bitangent and uv vectors 
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);

// Hit record struct ->
struct HitRecord
{
    float HitDistance;
    vec3 HitNormal;
};

// Gets normal for the intersected cube ->
vec3 GetNormal(vec3 v)
{
    vec3 s = sign(v);
    vec3 a = abs(v);
    vec3 n = mix(mix(vec3(0.0, 0.0, s.z), vec3(s.x, 0.0, 0.0), step(a.z, a.x)),
                 mix(vec3(0.0, s.y, 0.0), vec3(s.x, 0.0, 0.0), step(a.y, a.x)),
                 step(a.z, a.y));
    return n;
}

// Ray box intersection test
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

// Basic intersection test
HitRecord IntersectItemCube(vec3 r0, vec3 rD, vec3 BoxMinCoord, vec3 BoxMaxCoord)
{
    HitRecord result;
    result.HitDistance = RayBoxIntersectionTest(r0, rD, BoxMinCoord, BoxMaxCoord);
    float Transversal = step(result.HitDistance, INF);
    result.HitNormal = mix(-rD, GetNormal(r0 + rD * result.HitDistance - (BoxMinCoord + BoxMaxCoord) / 2.0f), Transversal);
    return result;
}

// Gets the percieved brightness of a pixel
float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}


// Remaps a value from one range to another
float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

// 3D rotation matrix 
// Todo : do this on the CPU
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

// RNG ->
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

// Samples ggx vndf ->
vec3 ImportanceSampleGGX(vec3 N, float roughness, vec2 Xi)
{
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
	
    float phi = 2.0 * PI * Xi.x;
    float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (alpha2 - 1.0) * Xi.y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	
    vec3 H;
    H.x = cos(phi) * sinTheta;
    H.y = sin(phi) * sinTheta;
    H.z = cosTheta;
	
    vec3 up = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(up, N));
    vec3 bitangent = cross(N, tangent);
	
    vec3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
    return normalize(sampleVec);
} 

// Gets a reflection direction ->
// (Picks a direction that is nearest to the normal after x samples)
vec3 GetReflectionDirection(vec3 N, float R) 
{
	R = max(R, 0.04f);
	float NearestDot = -100.0f;
	vec3 BestDirection;

	for (int i = 0 ; i < 2 ; i++) {
		vec2 Xi = hash2();
		Xi = Xi * vec2(1.0f, 0.75f);
		vec3 ImportanceSampled = ImportanceSampleGGX(N, R, Xi);
		float d = dot(ImportanceSampled,N);
		if (d > NearestDot) {
			BestDirection = ImportanceSampled;
			NearestDot = d;
		}
	}

	return BestDirection;
}

// Integrates ggx specular ->
vec3 SampleSpecular(vec3 I, vec3 N, float R, bool Aliased) {
    
    int Samples = int(mix(8, 32, clamp(R * 4.0f,0.,1.)));
    float x = mix(float(Samples) / max(float(u_AntialiasLevel) / 1.0f, 1.0f),float(Samples),float(!Aliased));
    Samples = int(floor(x));
    Samples = clamp(Samples, 2, 32);
    Samples += u_SpecularSampleBias;
    Samples = clamp(Samples, 2, 64);
    R = max(R, 0.05f);
    R = min(R, 0.5f);

    vec3 Specular = vec3(0.0f);

    for (int s = 0 ; s < Samples ; s++) {
        vec3 Microfacet = GetReflectionDirection(N, R);
        vec3 Reflected = reflect(I, Microfacet);
        Specular += texture(u_RandomHDRI, Reflected).xyz;
    }

    Specular /= float(Samples);

    return Specular;
}

// Fresnel schlick, with a roughness weight 
vec3 fresnelroughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
{
    float cosTheta = max(dot(norm, Eye), 0.0);
    const float magic = 2.4f;
    return F0 + (max(vec3(pow(1.0f - roughness, magic)), F0) - F0) * pow(clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
}

mat3 g_RotationMatrix;


// Intersects the item box
bool Intersect(vec2 LocalUV, vec2 ActualUV, out vec3 Normal) 
{
    // This shit took me so long to get right :(
    vec2 NDC = LocalUV * 2.0f - 1.0f;
    vec2 Clip = NDC;
    //float SinTime = sin(u_Time);
    
    vec3 RayDirection = g_RotationMatrix * vec3(NDC, 1.0f);
    vec3 RayOrigin = g_RotationMatrix * vec3(0.0f, 1.0f, -1.75f);

    RayDirection += vec3(0.0f, -0.52f, 0.0f);
    RayDirection = normalize(RayDirection);
    HitRecord result = IntersectItemCube(RayOrigin, RayDirection, vec3(0.0f), vec3(1.0f));

    if(result.HitDistance >= INF - 0.025f)
    {
        Normal = result.HitNormal;
        return false;
    }

    return true;
}

// basic tonemap used for testing 
vec3 BasicTonemap(vec3 color)
{
    float l = length(color);
    color = mix(color, color * 0.5f, l / (l + 1.0f));
    color = (color / sqrt(color * color + 1.0f));
    return color;
}

// Reinhard white shift preserving tonemap
float Reinhard2(float x)
{
    const float L_white = 4.0;
    return (x * (1.0 + x / (L_white * L_white))) / (1.0 + x);
}

// Reinhard white shift preserving tonemap
vec3 Reinhard2(vec3 x) 
{
    const float L_white = 4.0;
    return (x * (1.0 + x / (L_white * L_white))) / (1.0 + x);
}

// Reduces texture aliasing 
vec4 PixelArtFiltering(sampler2DArray tex, vec3 uvw, float LOD)
{
    vec2 uv = uvw.xy;
    vec2 res = vec2(textureSize(tex,0));
    uv = uv*res;
    vec2 seam = floor(uv+0.5);
    uv = seam + clamp( (uv-seam)/fwidth(uv), -0.5, 0.5);
    return textureLod(tex, vec3(uv/res, uvw.z), LOD);
}

// Saturation
vec3 BasicSaturation(vec3 Color, float Adjustment)
{
    const vec3 LuminosityCoefficients = vec3(0.2125f, 0.7154f, 0.0721f);
    vec3 Luminosity = vec3(dot(Color, LuminosityCoefficients));
    return mix(Luminosity, Color, Adjustment);
}

// Inverse fresnel 
float InverseSchlick(float f0, float VoH) 
{
    return 1.0 - clamp(f0 + (1.0f - f0) * pow(1.0f - VoH, 5.0f), 0.0f, 1.0f);
}

// Prevents bilinear artifacts by a bit 
vec4 TextureSmooth(sampler2DArray samp, vec3 uvw, float LOD) 
{
    vec2 uv = uvw.xy;
    vec2 textureResolution = textureSize(samp, 0).xy;
	uv = uv*textureResolution + 0.5f;
	vec2 iuv = floor(uv);
	vec2 fuv = fract(uv);
	uv = iuv + fuv*fuv*(3.0f-2.0f*fuv); 
	uv = (uv - 0.5f) / textureResolution;
	return textureLod(samp, vec3(uv,uvw.z), LOD).xyzw;
}

// 2D rotation matrix 
mat2 GenerateRotationMatrix2D(float angle)
{
	angle *= PI / 180.0;
    float s = sin(angle), c = cos(angle);
    return mat2(c, -s, s, c);
}

float g_Exposure;
vec3 g_SkyLightingUP;

// Calculates the radiance for a particular UV
vec3 Radiance(vec2 LocalUV, vec2 ActualUV, bool AliasSample) 
{
    // Convert to clip space 
    vec2 NDC = LocalUV * 2.0f - 1.0f;
    vec2 Clip = NDC;

    //float SinTime = sin(u_Time);

    // Generate matrices ->
    vec3 RayDirection = g_RotationMatrix * vec3(NDC, 1.0f);
    vec3 RayOrigin = g_RotationMatrix * vec3(0.0f, 1.0f, -1.75f);

    // Offset direction ->
    RayDirection += vec3(0.0f, -0.52f, 0.0f);

    RayDirection = normalize(RayDirection);

    // Intersect cube ->
    HitRecord result = IntersectItemCube(RayOrigin, RayDirection, vec3(0.0f), vec3(1.0f));

    // No hit ->
    if(result.HitDistance >= INF - 0.025f)
    {
        return imageLoad(o_OutputImage, ivec2(floor(ActualUV*u_Dimensions))).xyz;
    }

    // Hit occured, fuck me :) 
    vec3 IntersectionPoint = RayOrigin + RayDirection * result.HitDistance;
    vec3 IntersectionNormal = normalize(result.HitNormal);

    // Calc vectors ->
    vec2 UV;
    vec3 Tangent;
    vec3 Bitangent;
    CalculateVectors(IntersectionPoint, IntersectionNormal, Tangent, Bitangent, UV);
    UV.y = 1.0f - UV.y;

    // Generate TBN
    mat3 TBN = mat3(Tangent, Bitangent, IntersectionNormal);

    // Fetch ids ->
    ivec3 texture_ids;
    texture_ids.x = BlockAlbedoData[u_HeldBlockID];
    texture_ids.y = BlockNormalData[u_HeldBlockID];
    texture_ids.z = BlockPBRData[u_HeldBlockID];

    // Multitexturing hack ->
    if (u_HeldBlockID == u_GrassBlockProps[0])
	{
		if (IntersectionNormal == NORMAL_LEFT || IntersectionNormal == NORMAL_RIGHT || IntersectionNormal == NORMAL_FRONT || IntersectionNormal == NORMAL_BACK)
		{
			texture_ids.x = u_GrassBlockProps[4];
			texture_ids.y = u_GrassBlockProps[5];
			texture_ids.z = u_GrassBlockProps[6];
		}

		else if (IntersectionNormal == NORMAL_TOP)
		{
			texture_ids.x = u_GrassBlockProps[1];
			texture_ids.y = u_GrassBlockProps[2];
			texture_ids.z = u_GrassBlockProps[3];
		}

		else if (IntersectionNormal == NORMAL_BOTTOM)
		{
			texture_ids.x = u_GrassBlockProps[7];
			texture_ids.y = u_GrassBlockProps[8];
			texture_ids.z = u_GrassBlockProps[9];
		}
	}

    // fetch ->
    vec3 Albedo = PixelArtFiltering(u_AlbedoTextures, vec3(UV, float(texture_ids.x)), 3).xyz; // 512, 256, 128, 64
    vec3 Normal = PixelArtFiltering(u_NormalTextures, vec3(UV, float(texture_ids.y)), 3).xyz * 2.0f - 1.0f;
    Albedo = BasicSaturation(Albedo, 0.9f)*1.2f;
    vec4 PBR = PixelArtFiltering(u_PBRTextures, vec3(UV, float(texture_ids.z)), 3).xyzw;
    float IndirectOcclusion = pow(PBR.w,1.45f);
    
    if (PBR.y > 0.025f) {
        Albedo = BasicSaturation(Albedo, 0.8f);
    }

    int EmissivityMapID = BlockEmissiveData[u_HeldBlockID];
    bool BlockIsEmissive = EmissivityMapID > -0.5f;
    float Emissivity = BlockIsEmissive ? textureLod(u_EmissiveTextures, vec3(UV, float(EmissivityMapID)), 3).x : 0.0f;
    
    // Convert to tangent basis ->
    Normal = TBN * Normal;

    const vec3 VirtualLightPosition = normalize(-vec3(0.9f, -1.0f, 0.9f));
    //float HammonBRDF = DiffuseHammon(Normal, -RayDirection, VirtualLightPosition, PBR.x);
    //vec3 HammonDiffuse = HammonBRDF * Albedo;

    // Integrate lambert brdf ->
    float LambertianBRDF = clamp(dot(Normal, VirtualLightPosition) / PI, 0.0f, 1.0f);
    vec3 LambertianLighting = LambertianBRDF * LambertianBRDF * Albedo * vec3(1.1f, 1.1f, 1.1f);
    LambertianLighting = clamp(LambertianLighting, 0.0f, 1.0f);

    // If simple lighting is enabled, just use the lambert brdf
    if (u_SimpleLighting) {
        vec3 FakeAmbientShading = BasicSaturation(g_SkyLightingUP,0.75f)*0.4f;
        vec3 IntegratedBasic = (LambertianLighting+FakeAmbientShading*BasicSaturation(Albedo,1.2f))*g_Exposure*2.0f*max(Emissivity*4.0f,1.0f);
        return Reinhard2(IntegratedBasic*1.05f);
    }

    // Integrate direct lighting ->
    vec3 IndirectDiffuse = texture(u_RandomHDRIDiffuse, Normal).xyz;
    vec3 IndirectSkyModifier = texture(u_Sky, Normal).xyz;
    IndirectDiffuse = IndirectDiffuse * 0.75f;
    IndirectDiffuse *= 1.05f;
    IndirectDiffuse *= Albedo;

    vec3 I = RayDirection;

    vec3 IndirectSpecular = SampleSpecular(I, Normal, PBR.x, AliasSample);

    vec3 F0 = mix(vec3(0.04), Albedo, PBR.g);

    // Calculate color shift ->
    vec3 ColorBias = mix(g_SkyLightingUP*6.2f, vec3(0.9f), u_SunVisibility);
    ColorBias = mix(ColorBias,vec3(1.),BlockIsEmissive?(Emissivity>1.05f?0.75f:0.5f):0.0f);
    vec3 IndirectColorBias = mix(ColorBias, vec3(1.0f), 0.25f);
    vec3 DirectColorBias = ColorBias*1.5f;

    // Combine based on fresnel ->
    vec3 SpecularFactor = fresnelroughness(-I, Normal.xyz, vec3(F0), PBR.x) * 0.85f; 

    // Integrate indirect ->
    vec3 IntegratedIndirect = (IndirectDiffuse * (1.-SpecularFactor)) + (IndirectSpecular * SpecularFactor);
    IntegratedIndirect *= IndirectOcclusion;
    IntegratedIndirect *= IndirectColorBias;
    LambertianLighting *= DirectColorBias;

    // Combine direct and indirect ->
    vec3 Integrated = LambertianLighting + IntegratedIndirect;

    // Apply exposure and fake bloom

    vec3 AverageColorApproximate = TextureSmooth(u_AlbedoTextures, vec3(vec2(0.5f), float(texture_ids.x)), 8.0f).xyz+TextureSmooth(u_AlbedoTextures, vec3(vec2(0.8f), float(texture_ids.x)), 8.0f).xyz;
    AverageColorApproximate /= 2.0f;

    Emissivity = max(1.,Emissivity*3.5f);

    Integrated = Integrated * g_Exposure * Emissivity;
    Integrated *= 2.2f;

    // Fake bloom ->
    if (BlockIsEmissive) {
        float FakeGlow = TextureSmooth(u_EmissiveTextures, vec3(UV, float(EmissivityMapID)), 6.25f).x;
        vec3 FakeColor = pow(AverageColorApproximate, vec3(1.414f));
        Integrated += FakeGlow * (1.0f / g_Exposure) * mix(6.0f, 6.0f*6.0f, 1.-clamp01(u_SunVisibility)) * FakeColor;
    }

    // Tonemap and return

    const vec3 FakeBlueShift = vec3(0.925f,0.89f,1.0f);
    vec3 Tonemapped = Reinhard2(Integrated*FakeBlueShift);

    return Tonemapped ;
}

float Manhattan(vec2 p1, vec2 p2) {
    return abs(p1.x-p2.x)+abs(p1.y-p2.y);
}

// Entry ->
void main()
{
    if (!(v_TexCoords.x > 0.85f && v_TexCoords.y > 0.79f)) {
        return;
    }

    if (v_TexCoords.x > 0.972f || v_TexCoords.y > 0.972f) {
        return;
    }

    float DimensionRatio = u_Dimensions.x / u_Dimensions.y;

    // Compute globals ->

    // cube rotation matrix ->
    g_RotationMatrix = RotationMatrix(vec3(1.1f, 3.0f, 1.1f), u_Time * 0.75f);

    // Camera exposure ->
    g_Exposure = mix(4.0f, 12.0f, u_SunVisibility);

    // Basic sky ambient ->
    g_SkyLightingUP = texture(u_Sky,vec3(0.,1.,0.)).xyz;
    
    // RNG ->
    HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 200.0 * 20.0f;


    vec2 Texel = 1.0f / u_Dimensions;
    
    vec3 CenterRadiance;

    {
        vec2 TexCoords = v_TexCoords;
        vec2 LocalUV;
        LocalUV.x = remap(TexCoords.x, 0.75f, 1.0f, 0.0f, 1.0f);
        LocalUV.y = remap(TexCoords.y, 0.65f, 1.0f, 0.0f, 1.0f);
        CenterRadiance = Radiance(LocalUV, TexCoords, false).xyz;
    }


    // Custom Antialiasing Algorithm ->
    // -> Use screen space derivatives to find out luminance jumps
    // And only apply anti aliasing where it jumps. 
    // The anti aliasing algorithm samples multiple subpixels linearly from -0.5 -> +0.5

    float Threshold = mix(0.020f, 0.125f, float(u_SunVisibility));
    bool Aliased = fwidth(GetLuminance((CenterRadiance.xyz))) > Threshold;

    if (!Aliased || (v_TexCoords.x < 0.86f && v_TexCoords.y < 0.81f) || u_AntialiasLevel == 0) 
    {
        imageStore(o_OutputImage, ivec2(gl_FragCoord.xy), vec4(CenterRadiance, 1.0f));
        return;
    }

    #ifdef SMART_ANTI_ALIAS_DEBUG
        imageStore(o_OutputImage, ivec2(gl_FragCoord.xy), vec4(1.,0.,0.,1.));
        return;
    #endif

    // Pixel is aliased. Sample subpixels ->

    float AntiAliasingSamples = float(u_AntialiasLevel); // Need to make sure that the sample count is divisible by 2!
    AntiAliasingSamples = clamp(AntiAliasingSamples, 1, 16);
    float StepSize = 1.0f / AntiAliasingSamples;
    vec3 TotalRadiance = vec3(0.0f);
    
    for (float x = -0.5f; x < 0.5f; x += StepSize) 
    {
        for (float y = -0.5f; y < 0.5f; y += StepSize) 
        {
            if (Manhattan(vec2(x,y),vec2(0.0f)) < 0.001f) {
                TotalRadiance += CenterRadiance;
                continue;
            }

            vec2 TexCoords = v_TexCoords + vec2(x,y) * Texel;
            vec2 LocalUV;
            LocalUV.x = remap(TexCoords.x, 0.75f, 1.0f, 0.0f, 1.0f);
            LocalUV.y = remap(TexCoords.y, 0.65f, 1.0f, 0.0f, 1.0f);
            TotalRadiance += Radiance(LocalUV, TexCoords, true).xyz;
        }
    }
    
    TotalRadiance /= AntiAliasingSamples * AntiAliasingSamples;
    imageStore(o_OutputImage, ivec2(gl_FragCoord.xy), vec4(TotalRadiance,1.));
}


bool CompareVec3(vec3 v1, vec3 v2) {
	float e = 0.0125f;
	return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
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