// Generates textures gbuffer  


#version 430 core

#define EPS 0.000001f


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


layout (location = 0) out vec3 o_Albedo; // Albedo
layout (location = 1) out vec3 o_Normal; // High frequency normals
layout (location = 2) out vec4 o_PBR; // Roughness, metalness, displacement, emissivity
layout (location = 3) out float o_TextureAO; // Texture ambient occlusion

in vec2 v_TexCoords;

uniform sampler2D u_NonLinearDepth;
uniform sampler2D u_Normals;
uniform sampler2D u_BlockIDs;

uniform sampler2DArray u_BlockAlbedos;
uniform sampler2DArray u_BlockNormals;
uniform sampler2DArray u_BlockPBR;
uniform sampler2DArray u_BlockEmissive;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform int u_GrassBlockProps[10];

uniform int u_LavaBlockID;
uniform sampler3D u_LavaTextures[2];

uniform float u_Time;
uniform float uTime;
uniform bool u_UpdateGBufferThisFrame;

// POM 

uniform bool u_POM;
uniform bool u_HighQualityPOM;
uniform float u_POMHeight;

uniform mat4 u_View;


layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};


void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv, out vec4 UVDerivative);
vec4 GetTextureIDs(int BlockID, vec3 FlatNormal);


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

vec3 SampleNormalFromTex(sampler2D samp, vec2 txc) { 
    return GetNormalFromID(texture(samp, txc).x);
}

int GetBlockID(vec2 txc)
{
	float id = texture(u_BlockIDs, txc).r;
	return clamp(int(floor(id * 255.0f)), 0, 127);
}




vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).x;
	Dist = 1.0f / Dist;
	return vec4(u_InverseView[3].xyz + normalize(GetRayDirectionAt(txc)) * Dist, Dist);
}

bool CompareVec3(vec3 v1, vec3 v2) {
	float e = 0.0125f;
	return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
}


vec3 BasicTextureDistortion(vec3 UV) {
    vec2 UVxy = UV.xy;
    float time = u_Time;
    UV.x += sin(time * 0.25f);
    UV.y += pow(cos(time * 0.15f),2.);
    UV.x += cos(UV.x*10.0f + time)*0.3f;
    UV.y += sin(UV.y*5.0f + UV.x*4.0f + time*1.3f)*0.4f;
    UV.xy = mix(UV.xy,UVxy.xy,0.91f);
    return UV;
}


mat3 ArbitraryTBN(vec3 N) 
{
	vec3 up = vec3(0.,1.,0.);
	vec3 tangent = normalize(cross(up, N));
    vec3 bitangent = cross(N, tangent);
	return mat3(tangent, bitangent, N);
}


// Custom relief parallax mapping implementation 
// I think I fucked something up, it doesn't seem to work right 
// And there are some floating point artifacts 
// Todo : fix this motherfucker
vec3 ReliefParallax(vec3 WorldPosition, vec3 ViewVector, vec3 FlatNormal, mat3 BaseTangentBasis, vec2 FlatUV, int PBRID) {

	float ParallaxDepth = 0.2f;

	mat3 TangentBasis = BaseTangentBasis;//ArbitraryTBN(FlatNormal);

	// Convert view vector to tangent space ->
	vec3 TangentViewVector = TangentBasis * ViewVector;

	vec2 StartUV = clamp(FlatUV, EPS, 1.0f);
	vec2 TangentMarchDirection = TangentViewVector.xy; // Since it's in tangent space, all that matters is the x and y coordinate 
	TangentMarchDirection *= ParallaxDepth;
	TangentMarchDirection /= TangentViewVector.z;

	int MarchSteps = 64;
	int BinaryRefinementSteps = 8;

	float MarchStepSize = 1.0f / MarchSteps;
	float CurrentDepth = 1.0f;
	float BestDepth = 1.0f;



	// Find best depth ->
	for (int i = 1 ; i < MarchSteps ; i++) 
	{
		CurrentDepth -= MarchStepSize;
		
		vec2 SampleUV = StartUV + TangentMarchDirection * CurrentDepth;
		float SampleHeight = texture(u_BlockPBR, vec3(SampleUV, float(PBRID))).z;

		SampleHeight = 1.0f - SampleHeight;

		if (CurrentDepth >= SampleHeight) {
			BestDepth = CurrentDepth;
		}
	}

	CurrentDepth = BestDepth - MarchStepSize;


	// Perform a basic binary refinement step ->

	float BinaryStepSize = MarchStepSize;
	float CurrentBinaryDepth = CurrentDepth;
	float BestBinaryDepth = BestDepth;

	// Not really "step size" but I'm bad at naming variables :/

	const float BinaryRefineStepSize = 0.5f; 

	// Refinement loop
	for (int i = 1 ; i < BinaryRefinementSteps ; i++) {

		BinaryStepSize *= BinaryRefineStepSize;

		vec2 SampleUV = StartUV + TangentMarchDirection * CurrentDepth;
		float SampleHeight = texture(u_BlockPBR, vec3(SampleUV, float(PBRID))).z;

		SampleHeight = SampleHeight;

		if (CurrentBinaryDepth >= SampleHeight) {
			BestBinaryDepth = CurrentBinaryDepth;
			CurrentBinaryDepth -= BinaryStepSize * (1.0f / BinaryRefineStepSize);
		}

		CurrentBinaryDepth += BinaryStepSize;
	}

	// Calculate final coordinate 
	// Origin + Direction * BestDepth
	vec2 ParallaxCoordinate = StartUV + TangentMarchDirection * (BestBinaryDepth + 0.0001f);

	if (ParallaxCoordinate != clamp(ParallaxCoordinate, 0.00000001f, 1.0f)) {

		
	}

	return vec3(ParallaxCoordinate, BestBinaryDepth);


}


vec3 ParallaxGroundTruth(vec3 WorldPosition, vec3 ViewVector, vec3 FlatNormal, mat3 BaseTangentBasis, vec2 FlatUV, int PBRID) {
	

	mat3 TangentBasis = BaseTangentBasis;//ArbitraryTBN(FlatNormal);

	// Convert view vector to tangent space ->
	vec3 TangentViewVector = TangentBasis * ViewVector;

	vec2 StartUV = clamp(FlatUV, EPS, 1.0f);
	vec2 TangentMarchDirection = TangentViewVector.xy; // Since it's in tangent space, all that matters is the x and y coordinate 
	TangentMarchDirection *= 0.3f;
	TangentMarchDirection /= TangentViewVector.z;

	float BestHeight = 1.0f;

	float BaseHeight = texture(u_BlockPBR, vec3(FlatUV, float(PBRID))).z; 

	int StepCount = 1024;
	float StepSize = 1.0f / StepCount;

	vec2 CurrentSampleUV = FlatUV;

	float CurrentSampleHeight = BaseHeight;

	for (int Step = 0 ; Step < StepCount ; Step++) 
	{
	
		BestHeight -= StepSize; 
		CurrentSampleUV += TangentMarchDirection * StepSize; 
		CurrentSampleHeight = texture(u_BlockPBR, vec3(CurrentSampleUV, float(PBRID))).z; 

		if (BestHeight <= CurrentSampleHeight) {
			break;
		}
	}

	return vec3(CurrentSampleUV, BestHeight);

}

vec3 ParallaxOcclusion(vec3 WorldPosition, vec3 ViewVector, vec3 FlatNormal, mat3 BaseTangentBasis, vec2 FlatUV, int PBRID) {
	vec3 TangentViewVector = BaseTangentBasis * ViewVector;

	float NumLayers = 128; 
    float LayerDepth = 1.0 / NumLayers;

	float PomDepth = 0.2f;
  
	float CurrentLayerDepth = 0.0;
    
	vec2 P = TangentViewVector.xy * 1.0f; 
    
	vec2 DeltaTexCoords = P / NumLayers;
   
    vec2 InitialDeltaCoords = DeltaTexCoords;

    vec2  CurrentTexCoords = FlatUV;
    float CurrentDepthMapValue = texture(u_BlockPBR, vec3(FlatUV, float(PBRID))).z; CurrentDepthMapValue *= PomDepth; CurrentDepthMapValue = 1. - CurrentDepthMapValue;

	// Raymarch ->
    for (int i = 0 ; i < NumLayers ; i++)
    {
        if(CurrentLayerDepth < CurrentDepthMapValue)
        {
            CurrentTexCoords -= DeltaTexCoords;
            CurrentDepthMapValue = texture(u_BlockPBR, vec3(CurrentTexCoords, float(PBRID))).z; CurrentDepthMapValue *= PomDepth; CurrentDepthMapValue = 1. - CurrentDepthMapValue;
            CurrentLayerDepth += LayerDepth;
        }
    }

	// Interpolate ->
    vec2 PrevTexCoords = CurrentTexCoords + DeltaTexCoords;
    float AfterDepth  = CurrentDepthMapValue - CurrentLayerDepth;

	float d = texture(u_BlockPBR, vec3(PrevTexCoords, float(PBRID))).z ;

    float BeforeDepth = (1.0f-(d * PomDepth)) - CurrentLayerDepth + LayerDepth;

    float Weight = AfterDepth / (AfterDepth - BeforeDepth);
    vec2 FinalTexCoords = PrevTexCoords * Weight + CurrentTexCoords * (1.0 - Weight);

    return vec3(FinalTexCoords, CurrentLayerDepth);

}


vec3 Parallax(vec3 WorldPosition, vec3 ViewVector, vec3 FlatNormal, mat3 BaseTangentBasis, vec2 FlatUV, int PBRID) {
	
	if (!u_POM) {
		return vec3(FlatUV, EPS);
	}


	return ReliefParallax(WorldPosition, ViewVector, FlatNormal, BaseTangentBasis, FlatUV, PBRID);

}

void main() {
	vec3 PlayerPosition = u_InverseView[3].xyz;

	int BaseID = GetBlockID(v_TexCoords);

	bool ShouldUpdate = u_UpdateGBufferThisFrame || BaseID == u_LavaBlockID;

	if (!ShouldUpdate) {
		discard;
	}

	vec4 WorldPosition = GetPositionAt(u_NonLinearDepth, v_TexCoords);

	if (WorldPosition.w < 0.0f) {
		o_Albedo = vec3(0.0f);
		o_Normal = vec3(1.0f); 
		o_PBR = vec4(0.0f);
		o_TextureAO = 0.0f; 
		return;
	}

	vec3 FlatNormal = SampleNormalFromTex(u_Normals, v_TexCoords).rgb;
	vec4 data = GetTextureIDs(BaseID, FlatNormal);

	vec2 UV;
    vec3 Tangent, Bitangent;
    vec4 UVDerivative;
    CalculateVectors(WorldPosition.xyz, FlatNormal, Tangent, Bitangent, UV, UVDerivative); 

	vec2 tUV;
    vec3 tTangent, tBitangent;
    vec4 tUVDerivative;
    CalculateVectors(WorldPosition.xyz, abs(FlatNormal), tTangent, tBitangent, tUV, tUVDerivative); 


	bool IsLiquid =  (u_LavaBlockID == BaseID) ;
    bool IsLava =  (u_LavaBlockID == BaseID) ;
	mat3 tbn = mat3((Tangent), (Bitangent), (FlatNormal));
    vec3 DistortedUV = IsLiquid ? BasicTextureDistortion(vec3(UV,fract(u_Time*0.3f))) : vec3(UV,fract(u_Time*0.3f));
	
	vec3 ViewVector = normalize(WorldPosition.xyz - PlayerPosition);
	UV = IsLiquid ? DistortedUV.xy : Parallax(WorldPosition.xyz, ViewVector, FlatNormal, mat3(tTangent, tBitangent, abs(FlatNormal)), tUV.xy, int(data.z)).xy;

	UV = 1.0f - UV;

    vec3 NormalMapped = IsLava ? (texture(u_LavaTextures[1], DistortedUV).xyz) : (textureGrad(u_BlockNormals, vec3(UV, data.y), UVDerivative.xy, UVDerivative.zw).xyz);
    NormalMapped = NormalMapped * 2.0f - 1.0f;
    NormalMapped = tbn * NormalMapped;

    vec4 PBRMap = textureGrad(u_BlockPBR, vec3(UV, data.z), UVDerivative.xy, UVDerivative.zw).xyzw;
    float Emissivity = data.w > -0.5f ? texture(u_BlockEmissive, vec3(UV, data.w)).x : 0.0f;

    o_Normal = NormalMapped;
    o_PBR = vec4(PBRMap.x, PBRMap.y, PBRMap.z, Emissivity);

    o_PBR = clamp(o_PBR, 0.0f, 1.0f);

    o_TextureAO = PBRMap.w;

    o_TextureAO = clamp(o_TextureAO, 0.00000001f, 1.0f);

    vec3 AlbedoColor = IsLava ? (texture(u_LavaTextures[0], DistortedUV).xyz) : textureGrad(u_BlockAlbedos, vec3(UV, data.x), UVDerivative.xy, UVDerivative.zw).rgb;

    o_Albedo = AlbedoColor;

	const bool BloomLightLeakFix = true;
	if (BloomLightLeakFix&&!IsLiquid) {
		float lbiasx = 0.02f;
		float lbiasy = 0.02f;
		o_PBR.w *= float(UV.x > lbiasx && UV.x < 1.0f - lbiasx &&
							UV.y > lbiasy && UV.y < 1.0f - lbiasy);
	}
}

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

const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

vec4 GetTextureIDs(int BlockID, vec3 FlatNormal) 
{
	vec4 data = vec4(float(BlockAlbedoData[BlockID]),
				float(BlockNormalData[BlockID]),
				float(BlockPBRData[BlockID]),
				float(BlockEmissiveData[BlockID]));

    if (BlockID == u_GrassBlockProps[0])
	{
	    if (FlatNormal == NORMAL_LEFT || FlatNormal == NORMAL_RIGHT || FlatNormal == NORMAL_FRONT || FlatNormal == NORMAL_BACK)
	    {
	        data.x = u_GrassBlockProps[4];
	        data.y = u_GrassBlockProps[5];
	        data.z = u_GrassBlockProps[6];
	    }

	    else if (FlatNormal == NORMAL_TOP)
	    {
	        data.x = u_GrassBlockProps[1];
	        data.y = u_GrassBlockProps[2];
	        data.z = u_GrassBlockProps[3];
	    }

	    else if (FlatNormal == NORMAL_BOTTOM)
	    {
	        data.x = u_GrassBlockProps[7];
	        data.y = u_GrassBlockProps[8];
	        data.z = u_GrassBlockProps[9];
	    }
	}

	return data;
}