#version 430 core 
#define INF 100000.0f
#define PI 3.14159265359f
#define TAU (2.0f * PI)
#define clamp01(x) (clamp(x,0.,1.))

layout(rgba16f, binding = 1) uniform image2D o_OutputImage;

in vec2 v_TexCoords;

uniform float u_Time;
uniform vec2 u_Dimensions;
uniform vec3 u_SunDirection;
uniform vec3 u_MoonDirection;
uniform samplerCube u_Sky;
uniform samplerCube u_RandomHDRI;
uniform samplerCube u_RandomHDRIDiffuse;

// Block textures ->
uniform sampler2DArray u_AlbedoTextures;
uniform sampler2DArray u_NormalTextures;
uniform sampler2DArray u_EmissiveTextures;
uniform sampler2DArray u_PBRTextures;

uniform int u_HeldBlockID;

uniform bool u_Antialias;

uniform float u_SunVisibility;
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

const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);


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

HitRecord IntersectItemCube(vec3 r0, vec3 rD, vec3 BoxMinCoord, vec3 BoxMaxCoord)
{
    HitRecord result;
    result.HitDistance = RayBoxIntersectionTest(r0, rD, BoxMinCoord, BoxMaxCoord);

    float Transversal = step(result.HitDistance, INF);
    
    result.HitNormal = mix(-rD, GetNormal(r0 + rD * result.HitDistance - (BoxMinCoord + BoxMaxCoord) / 2.0f), Transversal);
    return result;
}


float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

vec3 TemperatureToRGB(float temperatureInKelvins);

vec3 SampleSunColor()
{
    const vec3 TemperatureModifier = TemperatureToRGB(5778.0f);
    vec3 SunTransmittance = texture(u_Sky, u_SunDirection.xyz).xyz;
    vec3 SunColor = SunTransmittance;
    SunColor *= TemperatureModifier;
    return SunColor * PI * 2.2f * 1.25f;
}

vec3 SampleMoonColor()
{
    vec3 MoonTransmittance = texture(u_Sky, u_MoonDirection).xyz;
    vec3 MoonColor = MoonTransmittance;
    MoonColor = MoonColor * PI * 1.25f;
    MoonColor = mix(MoonColor, vec3(GetLuminance(MoonColor)), 0.05f); 
    return MoonColor * 0.55f;
}

void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);


float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

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

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

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

vec3 SampleSpecular(vec3 I, vec3 N, float R) {
    
    int Samples = int(mix(12, 24, R * 2.0f));
    Samples = clamp(Samples, 2, 32);
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

vec3 fresnelroughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
{
    float cosTheta = max(dot(norm, Eye), 0.0);
    const float magic = 2.4f;
    return F0 + (max(vec3(pow(1.0f - roughness, magic)), F0) - F0) * pow(clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
}

bool Intersect(vec2 LocalUV, vec2 ActualUV, out vec3 Normal) 
{
    vec2 NDC = LocalUV * 2.0f - 1.0f;
    vec2 Clip = NDC;
    float SinTime = sin(u_Time);
    mat3 RotationMatrix = RotationMatrix(vec3(1.1f, 3.0f, 1.1f), u_Time * 0.35f);
    vec3 RayDirection = RotationMatrix * vec3(NDC, 1.0f);
    vec3 RayOrigin = RotationMatrix * vec3(0.0f, 1.0f, -1.75f);
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

vec3 BasicTonemap(vec3 color)
{
    float l = length(color);
    color = mix(color, color * 0.5f, l / (l + 1.0f));
    color = (color / sqrt(color * color + 1.0f));

    return color;
}

vec4 BetterTexture(sampler2DArray samp, vec3 wwuv, float LOD)
{
    vec2 textureResolution = textureSize(samp, 0).xy;
    vec2 uv = wwuv.xy;
    uv = uv * textureResolution + 0.5;
    vec2 iuv = floor(uv);
    vec2 fuv = fract(uv);
    uv = iuv + fuv * fuv * (3.0 - 2.0 * fuv);
    uv = (uv - 0.5) / textureResolution;
    return textureLod(samp, vec3(uv, wwuv.z), LOD).xyzw;
}

vec3 Radiance(vec2 LocalUV, vec2 ActualUV) 
{
    vec2 NDC = LocalUV * 2.0f - 1.0f;
    vec2 Clip = NDC;


    float SinTime = sin(u_Time);

    // Generate matrices ->
    mat3 RotationMatrix = RotationMatrix(vec3(1.1f, 3.0f, 1.1f), u_Time * 0.35f);
    vec3 RayDirection = RotationMatrix * vec3(NDC, 1.0f);
    vec3 RayOrigin = RotationMatrix * vec3(0.0f, 1.0f, -1.75f);

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

    vec2 UV;
    vec3 Tangent;
    vec3 Bitangent;

    CalculateVectors(IntersectionPoint, IntersectionNormal, Tangent, Bitangent, UV);

    UV.y = 1.0f - UV.y;

    mat3 TBN = mat3(Tangent, Bitangent, IntersectionNormal);

    ivec3 texture_ids;
    texture_ids.x = BlockAlbedoData[u_HeldBlockID];
    texture_ids.y = BlockNormalData[u_HeldBlockID];
    texture_ids.z = BlockPBRData[u_HeldBlockID];

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

    vec3 Albedo = BetterTexture(u_AlbedoTextures, vec3(UV, float(texture_ids.x)), 3).xyz;
    vec3 Normal = BetterTexture(u_NormalTextures, vec3(UV, float(texture_ids.y)), 3).xyz * 2.0f - 1.0f;
    vec3 PBR = BetterTexture(u_PBRTextures, vec3(UV, float(texture_ids.z)), 3).xyz;
    float Emissivity = BetterTexture(u_EmissiveTextures, vec3(UV, float(BlockEmissiveData[u_HeldBlockID])), 3).x;
    
    Normal = TBN * Normal;

    const vec3 VirtualLightPosition = normalize(vec3(0.9f, -1.0f, 0.9f));
    const vec3 VirtualLightRadiance = vec3(4.5f);
    vec3 SimpleLambertianBRDF = (clamp(dot(Normal, VirtualLightPosition), 0.0001f, 1.0f) / PI) * VirtualLightRadiance;
    SimpleLambertianBRDF = SimpleLambertianBRDF * Albedo;

    vec3 IndirectDiffuse = texture(u_RandomHDRIDiffuse, Normal).xyz;

    IndirectDiffuse *= Albedo;

    vec3 I = RayDirection;

    vec3 IndirectSpecular = SampleSpecular(I, Normal, PBR.x);

    vec3 F0 = mix(vec3(0.04), Albedo, PBR.g);

    vec3 SpecularFactor = fresnelroughness(-I, Normal.xyz, vec3(F0), PBR.x); 

    vec3 Integrated = (IndirectDiffuse * (1.-SpecularFactor)) + (IndirectSpecular * SpecularFactor * 0.75f);

    Emissivity = max(1.,Emissivity*6.5f);

    float Exposure = mix(0.8f, 4.20f, u_SunVisibility);

    vec3 ColorBias = mix(vec3(0.6f, 0.6f, 1.0f), vec3(0.9f), u_SunVisibility);

    Integrated = Integrated * Exposure * Emissivity * ColorBias;

    return BasicTonemap(Integrated);
}



void main()
{
    

    if (!(v_TexCoords.x > 0.75f && v_TexCoords.y > 0.65f)) {
        return;
    }

    HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 200.0 * 20.0f;

    bool AntiAlias = u_Antialias;
    vec2 Texel = 1.0f / u_Dimensions;
    
    if (AntiAlias) {
    
        float A = 4.0f, s = 1.0f / A, x, y;
        vec3 TotalRadiance = vec3(0.0f);
    
        for (x = -0.5f; x < 0.5f; x += s) 
        {
            for (y = -0.5f; y < 0.5f; y += s) 
            {
                 vec2 TexCoords = v_TexCoords + vec2(x,y) * Texel;
                 vec2 LocalUV;
                 LocalUV.x = remap(TexCoords.x, 0.75f, 1.0f, 0.0f, 1.0f);
                 LocalUV.y = remap(TexCoords.y, 0.65f, 1.0f, 0.0f, 1.0f);
    
                 TotalRadiance += Radiance(LocalUV, TexCoords).xyz;
            }
        }
    
        TotalRadiance /= A * A;

      
        imageStore(o_OutputImage, ivec2(gl_FragCoord.xy), vec4(TotalRadiance,1.));
    
    }
    
    else 
    {
        vec2 TexCoords = v_TexCoords;
        vec2 LocalUV;
        LocalUV.x = remap(TexCoords.x, 0.75f, 1.0f, 0.0f, 1.0f);
        LocalUV.y = remap(TexCoords.y, 0.65f, 1.0f, 0.0f, 1.0f);
        vec3 RadianceFinal = Radiance(LocalUV, TexCoords).xyz;


        imageStore(o_OutputImage, ivec2(gl_FragCoord.xy), vec4(RadianceFinal,1.));
    }
   
    //imageStore(o_OutputImage, ivec2(gl_FragCoord.xy), IsAtEdge ? vec4(1.,0.,0.,0.) : vec4(0.,1.,0.,0.));
    
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