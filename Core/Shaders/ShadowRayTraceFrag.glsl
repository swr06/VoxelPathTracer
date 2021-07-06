#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

layout (location = 0) out float o_Shadow;

in vec2 v_TexCoords;

uniform sampler3D u_VoxelData;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2DArray u_AlbedoTextures;
uniform sampler2D u_PrevShadowFBO;
uniform sampler3D u_DistanceFieldTexture;

uniform bool u_DoFullTrace;
uniform mat4 u_ShadowProjection;
uniform mat4 u_ShadowView;

uniform vec3 u_LightDirection;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

bool IsInVolume(in vec3 pos)
{
    if (pos.x < 0.0f || pos.y < 0.0f || pos.z < 0.0f || 
        pos.x > float(WORLD_SIZE_X - 1) || pos.y > float(WORLD_SIZE_Y - 1) || pos.z > float(WORLD_SIZE_Z - 1))
    {
        return false;    
    }   

    return true;
}

float GetVoxel(ivec3 loc)
{
    if (IsInVolume(loc))
    {
        return texelFetch(u_VoxelData, loc, 0).r;
    }
    
    return 0.0f;
}

float ToConservativeEuclidean(float Manhattan)
{
	return Manhattan == 1 ? 1 : Manhattan * 0.57735026918f;
}

float GetDistance(ivec3 loc)
{
    if (IsInVolume(loc))
    {
         return (texelFetch(u_DistanceFieldTexture, loc, 0).r);
    }
    
    return -1.0f;
}

float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, float blockType) 
{
	vec3 initial_origin = origin;
	const float epsilon = 0.01f;
	bool Intersection = false;

	int MinIdx = 0;
	ivec3 RaySign = ivec3(sign(direction));

	int itr = 0;

	for (itr = 0 ; itr < 350 ; itr++)
	{
		ivec3 Loc = ivec3(floor(origin));
		
		if (!IsInVolume(Loc))
		{
			Intersection = false;
			break;
		}

		float Dist = GetDistance(Loc) * 255.0f; 

		int Euclidean = int(floor(ToConservativeEuclidean(Dist)));

		if (Euclidean == 0)
		{
			break;
		}

		if (Euclidean == 1)
		{
			// Do the DDA algorithm for one voxel 

			ivec3 GridCoords = ivec3(origin);
			vec3 WithinVoxelCoords = origin - GridCoords;
			vec3 DistanceFactor = (((1 + RaySign) >> 1) - WithinVoxelCoords) * (1.0f / direction);

			MinIdx = DistanceFactor.x < DistanceFactor.y && RaySign.x != 0
				? (DistanceFactor.x < DistanceFactor.z || RaySign.z == 0 ? 0 : 2)
				: (DistanceFactor.y < DistanceFactor.z || RaySign.z == 0 ? 1 : 2);

			GridCoords[MinIdx] += RaySign[MinIdx];
			WithinVoxelCoords += direction * DistanceFactor[MinIdx];
			WithinVoxelCoords[MinIdx] = 1 - ((1 + RaySign) >> 1) [MinIdx]; // Bit shifts (on ints) to avoid division

			origin = GridCoords + WithinVoxelCoords;
			origin[MinIdx] += RaySign[MinIdx] * 0.0001f;

			Intersection = true;
		}

		else 
		{
			origin += int(Euclidean - 1) * direction;
		}
	}

	if (Intersection)
	{
		normal = vec3(0.0f);
		normal[MinIdx] = -RaySign[MinIdx];
		blockType = GetVoxel(ivec3(origin));
		return blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
	}

	return -1.0f;
}


vec2 ReprojectShadow (in vec3 pos)
{
	vec4 Projected = u_ShadowProjection * u_ShadowView * vec4(pos, 1.0f);
	Projected.xyz /= Projected.w;
	Projected.xy = Projected.xy * 0.5f + 0.5f;

	return Projected.xy;
}

void main()
{
	vec4 RayOrigin = texture(u_PositionTexture, v_TexCoords).rgba;
	vec3 RayDirection = normalize(u_LightDirection - (u_LightDirection * 0.1f));
	vec3 SampledNormal = texture(u_NormalTexture, v_TexCoords).rgb;
	vec3 Bias = SampledNormal * vec3(0.055f);
	
	if (u_DoFullTrace)
	{
		float T = -1.0f;

		float block_at = GetVoxel(ivec3(floor(RayOrigin.rgb + Bias)));
		 
		if (RayOrigin.w > 0.0f) 
		{
			float b; vec3 n;
			T = VoxelTraversalDF(RayOrigin.rgb + Bias, RayDirection, n, b);
		}

		if (T > 0.0f || block_at > 0) 
		{ 
			o_Shadow = 1.0f; 
		}
		
		else 
		{ 
			o_Shadow = 0.0f; 
		}
	}

	else 
	{
		vec2 PreviousFrameReprojected = ReprojectShadow(RayOrigin.xyz);

	    if (PreviousFrameReprojected.x > 0.0f && PreviousFrameReprojected.x < 1.0f && PreviousFrameReprojected.y > 0.0f && PreviousFrameReprojected.y < 1.0f)
		{
			o_Shadow = texture(u_PrevShadowFBO, PreviousFrameReprojected).r;
			return;
		}

		else 
		{
			float T = -1.0f;
			 
			if (RayOrigin.w > 0.0f) 
			{
				float b; vec3 n;
				T = VoxelTraversalDF(RayOrigin.rgb + Bias, RayDirection, n, b);
			}

			if (T > 0.0f) 
			{ 
				o_Shadow = 1.0f; 
			}
			
			else 
			{ 
				o_Shadow = 0.0f; 
			}
		}
	}
}