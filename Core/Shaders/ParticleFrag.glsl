// Shader for particles

#version 330 core

layout (location = 0) out vec4 o_Color;

in float v_Alpha;
in vec2 v_TexCoords;
in float v_Z;
in float v_IDX;
in vec3 v_BlockPosition;
in vec3 v_WorldPosition;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_ShadowTexture;
uniform sampler2D u_DiffuseTexture;
uniform sampler2D u_DiffuseTextureYoCoCg;
uniform sampler2DArray u_BlockTextures;

uniform vec2 u_Dimensions;
uniform vec3 u_SunDir;
uniform vec3 u_PlayerPos;
uniform mat4 u_CameraViewProjection;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

const vec3 SUN_COLOR_C = (vec3(192.0f, 216.0f, 255.0f) / 255.0f) * 6.5f;
const vec3 NIGHT_COLOR_C  = (vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 2.0f; 

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


bool GetPlayerIntersect(in vec3 WorldPos, in vec3 d)
{
	float x = 0.5 - 0.3f;
	vec3 VP = u_PlayerPos + vec3(-x, -x, +x);
    float t = RayCapsuleIntersection(WorldPos, d, VP, VP + vec3(0.0f, 1.0f, 0.0f), 0.5f);
    return t > 0.0f;
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

vec2 GetBlockSamplePosition()
{
	vec3 BlockPosition = v_BlockPosition;
//	BlockPosition -= vec3(0.5f);
	vec4 ProjectedPosition = u_CameraViewProjection * vec4(BlockPosition, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xy * 0.5f + 0.5f; // Convert to screen space!
}

vec3 GetPositionAt(vec2 ScreenSpace)
{
	vec2 Position = ScreenSpace * 2.0f - 1.0f;
	vec4 clip = vec4(Position.xy, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	vec3 Dir = vec3(u_InverseView * eye);
	vec3 Orig = u_InverseView[3].xyz;
	float Distance = texture(u_PositionTexture, ScreenSpace).r;
	Dir = normalize(Dir);
	return Orig + (Dir * Distance);
}

void main()
{
	vec3 SunDirection = normalize(u_SunDir);

	// Project ->
	vec2 ScreenSpaceCoordinates = gl_FragCoord.xy / u_Dimensions.xy;
	vec3 WorldPosition = GetPositionAt(ScreenSpaceCoordinates).xyz;
	vec4 ProjectedPosition = u_CameraViewProjection * vec4(WorldPosition, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;

	// Depth test
	if (ProjectedPosition.z < v_Z) 
	{
		discard;
	}

	vec3 MoonDirection = vec3(-SunDirection.x, -SunDirection.y, SunDirection.z);
	vec3 StrongerLightDirection = -SunDirection.y < 0.01f ? SunDirection : MoonDirection;

	float ParticleAlpha = 1.0f - v_Alpha;
	ParticleAlpha = exp(ParticleAlpha);

	vec3 AlbedoColor = texture(u_BlockTextures, vec3(v_TexCoords, v_IDX)).rgb; 
	float SunVisibility = clamp(dot(SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	
	// We need to approximate whether the block is in shadow and the gi from screen space info
	// So the following code is meant to find the best sample coordinate 
	float DistanceError = distance(WorldPosition, v_WorldPosition);
	const float ErrorThresh = sqrt(2.0f); // sqrt(2)
	vec2 SamplePosition = DistanceError < ErrorThresh ? ScreenSpaceCoordinates : GetBlockSamplePosition();
	
	// Intersect with player capsule
	float t1 = 0.0f, t2 = 0.0f;
    bool PlayerIntersect = GetPlayerIntersect(WorldPosition.xyz,SunDirection);

	float Shadow = PlayerIntersect ? 1.0f : clamp(texture(u_ShadowTexture, SamplePosition).r - 0.2f, 0.0000000000001f, 1.0f);
	vec4 DiffuseSample = texture(u_DiffuseTexture, SamplePosition).rgba;
	vec2 YoCoCg = texture(u_DiffuseTextureYoCoCg, SamplePosition).rg;
	
	// Sample spherical harmonic ->
	vec3 Diffuse = SHToIrridiance(DiffuseSample, YoCoCg);
	Diffuse *= 2.75f;
	
	// Fake direct 
	vec3 SUN_COLOR = Shadow > 0.01f ? vec3(0.0f) : SUN_COLOR_C;
	vec3 NIGHT_COLOR = Shadow > 0.01f ? vec3(0.0f) : NIGHT_COLOR_C;

	vec3 Color = vec3(1.0f);
	Color = (mix(SUN_COLOR, NIGHT_COLOR, SunVisibility) * AlbedoColor) + (AlbedoColor * Diffuse); 
	
	vec3 TonemappedColor = Color;
	
	// Tonemap ->
	
	//float l = dot(TonemappedColor, vec3(0.2126, 0.7152, 0.0722));
    //vec3 tc = TonemappedColor / (TonemappedColor + 1.0);
    //TonemappedColor = mix(TonemappedColor / (l + 1.0), tc, tc);
	
	TonemappedColor = 1.0f - exp(- Color);
	
	
	o_Color = vec4(TonemappedColor, clamp(ParticleAlpha, 0.0f, 1.0f));
}