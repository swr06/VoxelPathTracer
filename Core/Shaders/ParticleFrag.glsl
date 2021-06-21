#version 330 core

layout (location = 0) out vec4 o_Color;

in float v_Alpha;
in vec2 v_TexCoords;
in float v_Z;
in float v_IDX;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_ShadowTexture;
uniform sampler2D u_DiffuseTexture;
uniform sampler2DArray u_BlockTextures;

uniform vec2 u_Dimensions;
uniform vec3 u_SunDir;
uniform vec3 u_PlayerPos;
uniform mat4 u_CameraViewProjection;

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

void main()
{
	vec2 ScreenSpaceCoordinates = gl_FragCoord.xy / u_Dimensions.xy;
	vec3 WorldPosition = texture(u_PositionTexture, ScreenSpaceCoordinates).xyz;
	vec4 ProjectedPosition = u_CameraViewProjection * vec4(WorldPosition, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;

	if (ProjectedPosition.z < v_Z) 
	{
		discard;
	}

	vec3 MoonDirection = vec3(-u_SunDir.x, -u_SunDir.y, u_SunDir.z);
	vec3 StrongerLightDirection = -u_SunDir.y < 0.01f ? u_SunDir : MoonDirection;

	// Approximate lighting and gi using screen space coordinates (This is a huge approximation but looks fine for particles) 
	float ParticleAlpha = 1.0f - v_Alpha;
	ParticleAlpha = exp(ParticleAlpha * 0.9f);

	vec3 Color = texture(u_BlockTextures, vec3(v_TexCoords, v_IDX)).rgb; 
	float t1, t2;
    bool PlayerIntersect = RayBoxIntersect(u_PlayerPos + vec3(0.2f, 0.0f, 0.2f), u_PlayerPos - vec3(0.75f, 1.75f, 0.75f), WorldPosition, StrongerLightDirection, t1, t2);

	float SunVisibility = clamp(dot(u_SunDir, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	float Shadow = PlayerIntersect ? 1.0f : texture(u_ShadowTexture, ScreenSpaceCoordinates).r;
	vec3 Diffuse = texture(u_DiffuseTexture, ScreenSpaceCoordinates).rgb * 8.0f;
	
	vec3 SUN_COLOR = Shadow > 0.01f ? Diffuse : SUN_COLOR_C;
	vec3 NIGHT_COLOR = Shadow > 0.01f ? Diffuse : NIGHT_COLOR_C;
	Color = mix(Color * SUN_COLOR, Color * NIGHT_COLOR, SunVisibility); 
	o_Color = vec4(1.0f - exp(-Color), clamp(ParticleAlpha, 0.0f, 1.0f));
}