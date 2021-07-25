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
layout (location = 1) out float o_Normal;
layout (location = 2) out float o_BlockID;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform int u_CurrentFrame;

uniform sampler3D u_VoxelDataTexture;
uniform sampler3D u_DistanceFieldTexture;

uniform sampler2DArray u_AlbedoTextures;

uniform vec2 u_Dimensions;
uniform bool u_ShouldAlphaTest;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
};

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

vec3 MapSize = vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);
vec2 CalculateUV(vec3 world_pos, in vec3 normal);

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

int GetBlockID(float id)
{
	return clamp(int(floor(id * 255.0f)), 0, 127);
}

bool CompareVec3(vec3 v1, vec3 v2) {
	//float e = 0.0125f;
	//return abs(v1.x - v2.x) < e && abs(v1.y - v2.y) < e && abs(v1.z - v2.z) < e;
	return v1 == v2;
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


float GetNormalID(in vec3 normal)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );


	if (CompareVec3(normal, Normals[0]))
    {
        return 0.0f / 10.0f;
    }

    else if (CompareVec3(normal, Normals[1]))
    {
        return 1.0f / 10.0f;
    }

    else if (CompareVec3(normal, Normals[2]))
    {
        return 2.0f / 10.0f;
    }

    else if (CompareVec3(normal, Normals[3]))
    {
        return 3.0f / 10.0f;
    }
	
    else if (CompareVec3(normal, Normals[4]))
    {
        return 4.0f / 10.0f;
    }
    

    else if (CompareVec3(normal, Normals[5]))
    {
        return 5.0f / 10.0f;
    }

	return 0.0f;
}


float AlphaTest(float d, vec3 initial_origin, vec3 o, vec3 direction, inout vec3 normal, inout float blockType)
{
	vec2 SampleUV = CalculateUV(o, normal);
	SampleUV.y = 1.0f - SampleUV.y;
	int BlockID = GetBlockID(blockType);
	float Lod = textureQueryLod(u_AlbedoTextures, vec2(SampleUV)).x;
	float Alpha = textureLod(u_AlbedoTextures, vec3(SampleUV, float(BlockAlbedoData[BlockID])), clamp(Lod - 1.0f, 0.0f, 8.0f)).w;
	
	if (Alpha <= 0.02f) { 
		vec3 origin = o;
		const float epsilon = 0.01f;
		bool Intersection = false;

		ivec3 BaseGridCoords = ivec3(origin);

		int MinIdx = 0;
		ivec3 RaySign = ivec3(sign(direction));

		int itr = 0;

		for (itr = 0 ; itr < 50 ; itr++)
		{
			bool Continue = false;

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
				float BlockAt = GetVoxel(Loc);
				int IDAt = GetBlockID(BlockAt);

				if (IDAt != BlockID)
				{
					break;
				}

				Continue = true;
			}

			if (Euclidean == 1 || Continue)
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
			vec3 n = vec3(0.0f);
			n[MinIdx] = -RaySign[MinIdx];
			normal = (n);
			blockType = GetVoxel(ivec3(floor(origin)));
			return blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
		} 
	}

	return d;
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


	if (u_ShouldAlphaTest)
	{
		normal = vec3(0.0f);
		normal[MinIdx] = -RaySign[MinIdx];
		blockType = GetVoxel(ivec3(floor(origin)));
		float dd = blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
		return AlphaTest(dd, initial_origin.xyz, origin.xyz, direction, normal, blockType);
	}

	else 
	{
		if (Intersection)
		{
			normal = vec3(0.0f);
			normal[MinIdx] = -RaySign[MinIdx];
			blockType = GetVoxel(ivec3(floor(origin)));
			return blockType > 0.0f ? distance(origin, initial_origin) : -1.0f;
		}
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
		o_Normal = GetNormalID(normal);
	} 

	else
	{
		o_Normal = 0.0f;
	}

	o_HitDistance = t;

	bool transparent;

	o_BlockID = 0.0f;

	if (intersect)
	{
		o_BlockID = id;
	}
}

vec2 CalculateUV(vec3 world_pos, in vec3 normal)
{
	vec2 uv;

	const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
	const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
	const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
	const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
	const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
	const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

    if (CompareVec3(normal, NORMAL_TOP))
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (CompareVec3(normal, NORMAL_BOTTOM))
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (CompareVec3(normal, NORMAL_RIGHT))
    {
        uv = vec2(fract(world_pos.zy));
    }

    else if (CompareVec3(normal, NORMAL_LEFT))
    {
        uv = vec2(fract(world_pos.zy));
    }
    
    else if (CompareVec3(normal, NORMAL_FRONT))
    {
        uv = vec2(fract(world_pos.xy));
    }

     else if (CompareVec3(normal, NORMAL_BACK))
    {
        uv = vec2(fract(world_pos.xy));
    }

	return uv;
}