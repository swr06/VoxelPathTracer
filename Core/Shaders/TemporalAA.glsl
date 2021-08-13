// Possibly the worst temporal filter i've ever written 
// I'll need to change this shit soon


#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayOrigin;
in vec3 v_RayDirection;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_PositionTexture;

uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousPositionTexture;

uniform sampler2D u_CurrentNormalTexture;
uniform sampler2D u_PreviousNormalTexture;

uniform sampler2D u_CurrentBlockTexture;
uniform sampler2D u_PreviousBlockTexture;

uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

uniform bool u_Enabled;

vec2 View;
vec2 Dimensions;
vec2 TexCoord;

vec2 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

float FastDistance(in vec3 p1, in vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}

vec4 GetPositionAt(sampler2D pos_tex, vec2 txc)
{
	float Dist = texture(pos_tex, txc).r;
	return vec4(v_RayOrigin + normalize(v_RayDirection) * Dist, Dist);
}

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

bool Vec3Equal(vec3 a, vec3 b) {
	return (abs(a.x-b.x) < 0.01f) && (abs(a.y-b.y) < 0.01f) && (abs(a.z-b.z) < 0.01f);
}

vec3 NeighbourhoodClamping(vec3 color, vec3 tempColor) 
{
	vec2 neighbourhoodOffsets[8] = vec2[8]
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

	vec3 minclr = color, maxclr = color;

	for(int i = 0; i < 8; i++) 
	{
		vec2 offset = neighbourhoodOffsets[i] * View;
		vec3 clr = texture(u_CurrentColorTexture, TexCoord + offset, 0.0).rgb;
		minclr = min(minclr, clr);
		maxclr = max(maxclr, clr);

	}

	return clamp(tempColor, minclr - 0.00945f, maxclr + 0.00945f);
}



void main()
{
	Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
	View = 1.0f / Dimensions;
	TexCoord = v_TexCoords;

	vec3 CurrentColor = texture(u_CurrentColorTexture, TexCoord).rgb;

	if (!u_Enabled)
	{
		o_Color = CurrentColor;
		return;
	}

	vec4 WorldPosition = GetPositionAt(u_PositionTexture, v_TexCoords).rgba;

	if (WorldPosition.w < 0.0f)
	{
		o_Color = CurrentColor;
		return;
	}

	vec2 CurrentCoord = v_TexCoords;
	vec2 PreviousCoord = Reprojection(WorldPosition.xyz); 
	float bias = 0.01f;
	vec3 CurrentNormal = SampleNormalFromTex(u_CurrentNormalTexture, v_TexCoords).xyz;
	int CurrentBlock = clamp(int(floor((texture(u_CurrentBlockTexture, v_TexCoords).r) * 255.0f)), 0, 127);

	//PreviousCoord = texture(PreviousCoord, WorldPosition.xyz);

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		CurrentCoord.x > bias && CurrentCoord.x < 1.0f-bias &&
		CurrentCoord.y > bias && CurrentCoord.y < 1.0f-bias)
	{
		vec4 WorldPositionPrev = GetPositionAt(u_PreviousPositionTexture, PreviousCoord).rgba;
		vec3 PreviousNormal = SampleNormalFromTex(u_PreviousNormalTexture, PreviousCoord).xyz;
		float Error = distance(WorldPositionPrev.xyz, WorldPosition.xyz);
		vec3 PrevColor = texture(u_PreviousColorTexture, PreviousCoord).rgb;
		int PreviousBlock = clamp(int(floor((texture(u_PreviousBlockTexture, PreviousCoord).r) * 255.0f)), 0, 127);

		PrevColor = NeighbourhoodClamping(CurrentColor, PrevColor);

		if (WorldPositionPrev.w <= 0.0f || Error > 0.5f || !Vec3Equal(PreviousNormal, CurrentNormal) ||
			PreviousBlock != CurrentBlock)
		{
			o_Color = CurrentColor;
			return;
		}

		// Construct our motion vector
		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;
		float BlendFactor = exp(-length(velocity)) * 0.2f + 0.6f;
		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, clamp(BlendFactor, 0.025f, 0.925f));
	}

	else 
	{
		o_Color = CurrentColor;
	}
}

