#version 430 core

/*
Traversal Paper used : https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
*/

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

#define MULTIPLE_TEXTURING_GRASS
#define ALPHA_TESTING

layout (location = 0) out float o_HitDistance;
layout (location = 1) out vec3 o_Normal;
layout (location = 2) out vec4 o_Data;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform int u_CurrentFrame;

uniform sampler3D u_VoxelDataTexture;
uniform sampler3D u_DistanceFieldTexture;

uniform sampler2DArray u_AlbedoTextures;

uniform vec2 u_Dimensions;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

// Temporary solution to have multi texturing for grass blocks
// Data stored : 
// Grass block ID
// top face index (albedo, normal, pbr),
// right/left/front/back face index (albedo, normal, pbr),
// bottom face index (albedo, normal, pbr)
// 9 + 1 ints total

uniform int u_GrassBlockProps[10];

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

vec3 MapSize = vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);

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
        return texelFetch(u_VoxelDataTexture, loc, 0).r;
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

bool VoxelExists(in vec3 loc)
{
    if (GetVoxel(ivec3(loc)) > 0.0f) 
    {
        return true;
    }

    return false;
}

float VoxelTraversalDF(vec3 origin, vec3 direction, inout vec3 normal, inout float blockType) 
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
		blockType = GetVoxel(ivec3(floor(origin)));
		return blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
	}

	return -1.0f;
}



void main()
{
    Ray r;
    r.Origin = v_RayOrigin;
    r.Direction = normalize(v_RayDirection);

	int normal_idx = 0;

	vec3 normal;
	float id;


	float t = VoxelTraversalDF(r.Origin, r.Direction, normal, id);
	bool intersect = t > 0.0f && id > 0;

	if (intersect)
	{
		o_Normal = normal;
	} 

	else
	{
		o_Normal = vec3(-1.0f);
	}

	o_HitDistance = t;

	int reference_id;
	vec4 texture_ids;
	bool transparent;

	if (intersect)
	{
		reference_id = clamp(int(floor(id * 255.0f)), 0, 127);
		texture_ids.xyz = vec3(
			float(BlockAlbedoData[reference_id]),
			float(BlockNormalData[reference_id]),
			float(BlockPBRData[reference_id])
		);
		texture_ids.w = float(BlockEmissiveData[reference_id]);
	}

	else 
	{
		texture_ids = vec4(-1.0f);
	}

	const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
	const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
	const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
	const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
	const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
	const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

	#ifdef MULTIPLE_TEXTURING_GRASS

	// /* Specific to grass texture. Temporary solution */
	if (reference_id == u_GrassBlockProps[0])
	{
	    if (normal == NORMAL_LEFT || normal == NORMAL_RIGHT || normal == NORMAL_FRONT || normal == NORMAL_BACK)
		{
			texture_ids.x = u_GrassBlockProps[4];
			texture_ids.y = u_GrassBlockProps[5];
			texture_ids.z = u_GrassBlockProps[6];
		}

		else if (normal == NORMAL_TOP)
		{
			texture_ids.x = u_GrassBlockProps[1];
			texture_ids.y = u_GrassBlockProps[2];
			texture_ids.z = u_GrassBlockProps[3];
		}

		else if (normal == NORMAL_BOTTOM)
		{
			texture_ids.x = u_GrassBlockProps[7];
			texture_ids.y = u_GrassBlockProps[8];
			texture_ids.z = u_GrassBlockProps[9];
		}
	}

	#endif

	o_Data = vec4(texture_ids);
}