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

uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

uniform bool u_Enabled;
uniform bool u_BlockModified;

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
vec3 clipAABB(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

vec3 ClampColor(vec3 Color) 
{
    vec3 MinColor = vec3(100.0);
	vec3 MaxColor = vec3(-100.0); 
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture,0);

    for(int x = -2; x <= 2; x++) 
	{
        for(int y = -2; y <= 2; y++) 
		{
            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize).rgb; 
			float DepthAt = texture(u_PositionTexture, v_TexCoords).x;

			if (DepthAt <= 0.0f) { continue; }

            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

    return clipAABB(Color, MinColor, MaxColor);
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

	if (WorldPosition.w <= 0.0f)
	{
		o_Color = CurrentColor;
		return;
	}

	vec2 CurrentCoord = v_TexCoords;
	vec2 PreviousCoord = Reprojection(WorldPosition.xyz); 
	float bias = 0.01f;

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		CurrentCoord.x > bias && CurrentCoord.x < 1.0f-bias &&
		CurrentCoord.y > bias && CurrentCoord.y < 1.0f-bias)
	{
		vec3 PrevColor = texture(u_PreviousColorTexture, PreviousCoord).rgb;
		PrevColor = ClampColor(PrevColor);

		// Construct our motion vector
		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;
		float BlendFactor = exp(-length(velocity)) * 0.7f + 0.325f;
		BlendFactor = u_BlockModified ? 0.4f : BlendFactor;
		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, clamp(BlendFactor, 0.001f, 0.95f));
	}

	else 
	{
		o_Color = CurrentColor;
	}
}

