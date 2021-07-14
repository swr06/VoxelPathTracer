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

bool RayBoxIntersect(const vec3 boxMin, const vec3 boxMax, vec3 r0, vec3 rD, out float t_min, out float t_max) 
{
	vec3 inv_dir = 1.0f / rD;
	vec3 tbot = inv_dir * (boxMin - r0);
	vec3 ttop = inv_dir * (boxMax - r0);
	vec3 tmin = min(ttop, tbot);
	vec3 tmax = max(ttop, tbot);
	vec2 t = max(tmin.xx, tmin.yz);
	float t0 = max(t.x, t.y);
	t = min(tmax.xx, tmax.yz);
	float t1 = min(t.x, t.y);
	t_min = t0;
	t_max = t1;
	return t1 > max(t0, 0.0);
}

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

vec2 GetBlockSamplePosition()
{
	vec4 ProjectedPosition = u_CameraViewProjection * vec4(v_BlockPosition, 1.0f);
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
	vec2 ScreenSpaceCoordinates = gl_FragCoord.xy / u_Dimensions.xy;
	vec3 WorldPosition = GetPositionAt(ScreenSpaceCoordinates).xyz;
	vec4 ProjectedPosition = u_CameraViewProjection * vec4(WorldPosition, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;

	// Depth test
	if (ProjectedPosition.z < v_Z) 
	{
		discard;
	}

	vec3 MoonDirection = vec3(-u_SunDir.x, -u_SunDir.y, u_SunDir.z);
	vec3 StrongerLightDirection = -u_SunDir.y < 0.01f ? u_SunDir : MoonDirection;

	float ParticleAlpha = 1.0f - v_Alpha;
	ParticleAlpha = exp(ParticleAlpha);

	vec3 Color = texture(u_BlockTextures, vec3(v_TexCoords, v_IDX)).rgb; 
	float SunVisibility = clamp(dot(u_SunDir, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	
	// We need to approximate whether the block is in shadow and the gi from screen space info
	// So the following code is meant to find the best sample coordinate 
	float DistanceError = distance(WorldPosition, v_WorldPosition);
	const float ErrorThresh = floor(1.41421); // sqrt(2)
	vec2 SamplePosition = DistanceError < ErrorThresh ? ScreenSpaceCoordinates : GetBlockSamplePosition();
	
	// Apply the player shadow as well :p
	float t1 = 0.0f, t2 = 0.0f;
    bool PlayerIntersect = RayBoxIntersect(u_PlayerPos + vec3(0.2f, 0.0f, 0.2f), u_PlayerPos - vec3(0.75f, 1.75f, 0.75f), WorldPosition, StrongerLightDirection, t1, t2);

	float Shadow = PlayerIntersect ? 1.0f : texture(u_ShadowTexture, SamplePosition).r;
	vec4 DiffuseSample = texture(u_DiffuseTexture, SamplePosition).rgba;
	vec2 YoCoCg = texture(u_DiffuseTextureYoCoCg, SamplePosition).rg;
	vec3 Diffuse = SHToIrridiance(DiffuseSample, YoCoCg, vec3(0.0f, 1.0f, 0.0f));
	Diffuse *= 7.0f;
	vec3 SUN_COLOR = Shadow > 0.01f ? Diffuse : SUN_COLOR_C;
	vec3 NIGHT_COLOR = Shadow > 0.01f ? Diffuse : NIGHT_COLOR_C;

	Color = mix(Color * SUN_COLOR, Color * NIGHT_COLOR, SunVisibility); 
	o_Color = vec4(1.0f - exp(-Color), clamp(ParticleAlpha, 0.0f, 1.0f));
}