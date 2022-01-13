// Temporal Anti Aliasing implementation 
// Resources :
//	   TAA (Neighbourhood clamping, box clipping) : INSIDE TAA
//     Variance Clipping : https://www.dropbox.com/sh/dmye840y307lbpx/AAAQpC0MxMbuOsjm6XmTPgFJa	

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
uniform mat4 u_InversePrevProjection;
uniform mat4 u_InversePrevView;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_View;
uniform mat4 u_Projection;

uniform bool u_Enabled;
uniform bool u_BlockModified;

uniform vec3 u_CameraHistory[2];

vec2 View;
vec2 Dimensions;
vec2 TexCoord;

// reprojection.
vec2 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

// manhattan
float FastDistance(in vec3 p1, in vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
}


// current frame
vec3 GetRayDirection(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec4 SampleBasePosition(sampler2D pos_tex)
{
	float Dist = 1./texture(pos_tex, v_TexCoords).r;
	return vec4(v_RayOrigin + normalize(v_RayDirection) * Dist, Dist);
}



// prev position 
vec3 GetPREVRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InversePrevProjection * clip), -1.0, 0.0);
	return vec3(u_InversePrevView * eye);
}

vec4 GetPREVPositionAt(vec2 txc)
{
	float Dist = 1./texture(u_PreviousPositionTexture, txc).r;
	return vec4(u_InversePrevView[3].xyz + normalize(GetPREVRayDirectionAt(txc)) * Dist, Dist);
}




// AABB clipping - from inside.
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


vec3 rgb2ycocg(in vec3 rgb)
{
    float co = rgb.r - rgb.b;
    float t = rgb.b + co / 2.0;
    float cg = rgb.g - t;
    float y = t + cg / 2.0;
    return vec3(y, co, cg);
}


vec3 ycocg2rgb(in vec3 ycocg)
{
    float t = ycocg.r - ycocg.b / 2.0;
    float g = ycocg.b + t;
    float b = t - ycocg.g / 2.0;
    float r = ycocg.g + b;
    return vec3(r, g, b);
}

float Luminance(vec3 RGB )
{
    return dot(RGB, vec3(0.2126f, 0.7152f, 0.0722f));
}

vec3 Reinhard(vec3 RGB )
{
	RGB *= 2.44f;
    return vec3(RGB) / (vec3(1.0f) + Luminance(RGB));
}

vec3 InverseReinhard(vec3 RGB)
{
    return (RGB / (vec3(1.0f) - Luminance(RGB))) / 2.44f;
}


vec3 SampleHistory(vec2 Reprojected, vec4 WorldPosition) 
{
    vec3 MinColor = vec3(100.0);
	vec3 MaxColor = vec3(-100.0); 
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture,0);
	vec2 BestOffset = vec2(0.0f);
	float BestDistance = 1000.0f;
	int KernelX = 2;
	int KernelY = 2;
	bool FindBestPixel = false;

    for(int x = -KernelX; x <= KernelX; x++) 
	{
        for(int y = -KernelY; y <= KernelY; y++) 
		{
			if (WorldPosition.w > 0.0f&&FindBestPixel) {
				vec4 PrevPositionAt = GetPREVPositionAt(Reprojected + vec2(x, y) * TexelSize).xyzw;
				
				if (PrevPositionAt.w > 0.0f) 
				{
					float DistanceAt = FastDistance(WorldPosition.xyz, PrevPositionAt.xyz);
					if (DistanceAt < BestDistance) {
						BestDistance = DistanceAt;
						BestOffset = vec2(x,y);
					}
				}
			}


            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize).rgb; 

            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

	BestOffset = FindBestPixel ? BestOffset : vec2(0.0f);

	float Bias = 0.00001f;
	MinColor -= Bias;
	MaxColor += Bias;


	vec3 Color = texture(u_PreviousColorTexture, Reprojected ).xyz;
	
    return (clipAABB((Color), (MinColor), (MaxColor))).xyz;
    return (Color).xyz;
}



///
vec3 DiagonalColorClamp(vec2 Reprojected) {
	vec3 MinBox = vec3(100.0f);
	vec3 MaxBox = vec3(-100.0f);

	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);

	vec3 CenterColor = texture(u_CurrentColorTexture, Reprojected).xyz;

	const vec2 SampleOffsets[8] = vec2[8](vec2(-1, 1), vec2(0, 1), vec2(1, 1), vec2(-1, 0), vec2(1, 0), vec2(-1,-1), vec2(0,-1), vec2(1,-1));
	vec3 Samples[8];

	for (int x = 0 ; x < 8 ; x++) {
		vec2 Offset = SampleOffsets[x];
		vec2 Coord = Reprojected+Offset*TexelSize;
		vec3 SampleColor = texture(u_CurrentColorTexture, Coord).xyz; 
		MinBox = min(SampleColor, MinBox);
		MaxBox = min(SampleColor, MaxBox);
		Samples[x] = SampleColor;
	}

	// Diagonal Clamp
	const vec3 dgx = normalize(vec3(1.0f)); 
    const vec3 dgz = normalize(cross(dgx, vec3(0,1,0)));
    const vec3 dgy = cross(dgz, dgx);
    mat3x3 rgb2drgb = mat3x3(dgx, dgy, dgz);

    vec3 dc = rgb2drgb * CenterColor;
    vec3 mnd = dc;
	vec3 mxd = dc;

    for (int i = 0; i < 8; i++)
    {
        vec3 ds = rgb2drgb * Samples[i];
        mnd = min(mnd, ds);
        mxd = max(mxd, ds);
    }


	// History
	vec3 History = texture(u_PreviousColorTexture, Reprojected).xyz;
	History = clamp(History, MinBox, MaxBox);
	vec3 DiagonalTransform = rgb2drgb * History;
    DiagonalTransform = clamp(DiagonalTransform, mnd, mxd );
    return DiagonalTransform * rgb2drgb;
}
///

// inside taa 
vec3 clipToAABB(in vec3 cOld, in vec3 cNew, in vec3 centre, in vec3 halfSize)
{
    if (all(lessThanEqual(abs(cOld - centre), halfSize))) {
        return cOld;
    }
    
    vec3 dir = (cNew - cOld);
    vec3 near = centre - sign(dir) * halfSize;
    vec3 tAll = (near - cOld) / dir;
    float t = 1e20;
    for (int i = 0; i < 3; i++) {
        if (tAll[i] >= 0.0 && tAll[i] < t) {
            t = tAll[i];
        }
    }
    
    if (t >= 1e20) {
		return cOld;
    }
    return cOld + dir * t;
}

// Variance clipping
vec3 VarianceClip(vec2 Reproj, vec3 CenterCurrent) 
{
	vec2 Offsets[4] = vec2[4](vec2(-1.0,  0.0), 
								vec2( 1.0,  0.0), 
								vec2( 0.0, -1.0), 
								vec2( 0.0,  1.0));
	
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture, 0);
	vec3 C = texture(u_PreviousColorTexture, Reproj).xyz;
	vec3 Averaged = rgb2ycocg(CenterCurrent.rgb);
    vec3 StandardDeviation = Averaged * Averaged;

    for (int i = 0; i < 4; i++) 
	{
        vec3 Sample = rgb2ycocg(texture(u_CurrentColorTexture, (Reproj + Offsets[i] * TexelSize)).rgb);
        Averaged += Sample;
        StandardDeviation += Sample * Sample;
    }
	
    Averaged /= 5.0f;
    StandardDeviation = sqrt(StandardDeviation / 5.0f - Averaged * Averaged);
	vec3 ClippedColor = C;
	ClippedColor = ycocg2rgb(clipToAABB(rgb2ycocg(ClippedColor), rgb2ycocg(CenterCurrent.rgb), Averaged, StandardDeviation));
	return ClippedColor;
}

vec3 ToViewSpace(vec3 x) {
	return vec3(u_View * vec4(x, 1.0f));
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

	vec4 WorldPosition = SampleBasePosition(u_PositionTexture).rgba;

	
	vec3 rD = GetRayDirection(v_TexCoords).xyz;
	rD = normalize(rD);

	bool MotionVector = false;
	float ReduceWeight = 1.0f;

	if (WorldPosition.w < 0.001f)
	{
		const float rpd = 32.0f;
		WorldPosition.xyz = u_InverseView[3].xyz + rD * rpd;
		WorldPosition.w = rpd;
		MotionVector = true;
		ReduceWeight = 0.45f; // Accumulate fewer frames

	}


	vec2 CurrentCoord = v_TexCoords;
	vec2 PreviousCoord;
	
	vec3 PreviousCameraPos = u_CameraHistory[1];
	vec3 CurrentCameraPos = u_CameraHistory[0];
	vec3 Offset = (PreviousCameraPos - CurrentCameraPos);

	PreviousCoord = Reprojection(WorldPosition.xyz);
	
	
	//if (MotionVector) {
	//	
	//	vec2 ReprojectedA = Reprojection(PreviousCameraPos);
	//	vec2 ReprojectedB = Reprojection(CurrentCameraPos);
	//	PreviousCoord -= ReprojectedA - ReprojectedB;
	//}

	float bias = 0.01f;

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		CurrentCoord.x > bias && CurrentCoord.x < 1.0f-bias &&
		CurrentCoord.y > bias && CurrentCoord.y < 1.0f-bias)
	{
		// used for testing
		////vec3 PrevColor = AdvancedClamp(PreviousCoord); //SampleHistory(PreviousCoord, WorldPosition.xyzw);////
		////vec3 PrevColor = VarianceClip(PreviousCoord, CurrentColor.xyz);////

		vec3 PrevColor = SampleHistory(PreviousCoord, WorldPosition.xyzw);

		// Construct our motion vector
		vec2 velocity = (TexCoord - PreviousCoord.xy) * Dimensions;
		float BlendFactor = exp(-length(velocity)) * 0.7f + 0.45f;
		BlendFactor = u_BlockModified ? 0.1f : BlendFactor;
		o_Color = InverseReinhard(mix(Reinhard(CurrentColor.xyz), Reinhard(PrevColor.xyz), clamp(BlendFactor * ReduceWeight, 0.001f, 0.95f)));
	}

	else 
	{
		o_Color = CurrentColor;
	}
}

