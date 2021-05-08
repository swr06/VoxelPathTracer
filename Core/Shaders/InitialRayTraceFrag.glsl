#version 330 core

/*
Traversal Paper used : https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
*/

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

#define MULTIPLE_TEXTURING_GRASS

layout (location = 0) out vec4 o_Position;
layout (location = 1) out vec3 o_Normal;
layout (location = 2) out vec4 o_Data;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform int u_CurrentFrame;

uniform sampler3D u_VoxelDataTexture;
uniform vec2 u_Dimensions;
uniform vec3 BLOCK_TEXTURE_DATA[128];

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

const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);


float sum(vec3 v)
{
    return v.x + v.y + v.z;
}

bool IsInVoxelizationVolume(in vec3 pos)
{
    if (pos.x < 0.0f || pos.y < 0.0f || pos.z < 0.0f || 
        pos.x > float(WORLD_SIZE_X) || pos.y > float(WORLD_SIZE_Y) || pos.z > float(WORLD_SIZE_Z))
    {
        return false;    
    }   

    return true;
}

float GetVoxel(ivec3 loc)
{
    if (IsInVoxelizationVolume(loc))
    {
         return texelFetch(u_VoxelDataTexture, loc, 0).r;
    }
    
    return 0.0f;
}

bool VoxelExists(in vec3 loc)
{
    if (GetVoxel(ivec3(loc)) > 0.0f) 
    {
        return true;
    }

    return false;
}

// Calculates uv from world position
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv, out int NormalIndex)
{
    if (normal == NORMAL_TOP)
    {
        uv = vec2(fract(world_pos.xz));
		NormalIndex = 0;
    }

    else if (normal == NORMAL_BOTTOM)
    {
        uv = vec2(fract(world_pos.xz));
		NormalIndex = 1;
    }

    else if (normal == NORMAL_RIGHT)
    {
        uv = vec2(fract(world_pos.zy));
		NormalIndex = 2;
    }

    else if (normal == NORMAL_LEFT)
    {
        uv = vec2(fract(world_pos.zy));
		NormalIndex = 3;
    }
    
    else if (normal == NORMAL_FRONT)
    {
        uv = vec2(fract(world_pos.xy));
		NormalIndex = 4;
    }

     else if (normal == NORMAL_BACK)
    {
        uv = vec2(fract(world_pos.xy));
		NormalIndex = 5;
    }
}

float ProjectToCube(vec3 ro, vec3 rd) 
{	
	float tx1 = (0 - ro.x) / rd.x;
	float tx2 = (MapSize.x - ro.x) / rd.x;

	float ty1 = (0 - ro.y) / rd.y;
	float ty2 = (MapSize.y - ro.y) / rd.y;

	float tz1 = (0 - ro.z) / rd.z;
	float tz2 = (MapSize.z - ro.z) / rd.z;

	float tx = max(min(tx1, tx2), 0);
	float ty = max(min(ty1, ty2), 0);
	float tz = max(min(tz1, tz2), 0);

	float t = max(tx, max(ty, tz));
	
	return t;
}

float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType) 
{
	vec3 origin = orig;
	const float epsilon = 0.001f;
	float t1 = max(ProjectToCube(origin, direction) - epsilon, 0.0f);
	origin += t1 * direction;

	int mapX = int(floor(origin.x));
	int mapY = int(floor(origin.y));
	int mapZ = int(floor(origin.z));

	float sideDistX;
	float sideDistY;
	float sideDistZ;

	float deltaDX = abs(1.0f / direction.x);
	float deltaDY = abs(1.0f / direction.y);
	float deltaDZ = abs(1.0f / direction.z);
	float T = -1.0;

	int stepX;
	int stepY;
	int stepZ;

	int hit = 0;
	int side;

	if (direction.x < 0)
	{
		stepX = -1;
		sideDistX = (origin.x - mapX) * deltaDX;
	} 
	
	else 
	{
		stepX = 1;
		sideDistX = (mapX + 1.0 - origin.x) * deltaDX;
	}

	if (direction.y < 0) 
	{
		stepY = -1;
		sideDistY = (origin.y - mapY) * deltaDY;
	} 
	
	else 
	{
		stepY = 1;
		sideDistY = (mapY + 1.0 - origin.y) * deltaDY;
	}

	if (direction.z < 0) 
	{
		stepZ = -1;
		sideDistZ = (origin.z - mapZ) * deltaDZ;
	} 
	
	else 
	{
		stepZ = 1;
		sideDistZ = (mapZ + 1.0 - origin.z) * deltaDZ;
	}

	for (int i = 0; i < 175; i++) 
	{
		if ((mapX >= MapSize.x && stepX > 0) || (mapY >= MapSize.y && stepY > 0) || (mapZ >= MapSize.z && stepZ > 0)) break;
		if ((mapX < 0 && stepX < 0) || (mapY < 0 && stepY < 0) || (mapZ < 0 && stepZ < 0)) break;

		if (sideDistX < sideDistY && sideDistX < sideDistZ) 
		{
			sideDistX += deltaDX;
			mapX += stepX;
			side = 0;
		} 
		
		else if (sideDistY < sideDistX && sideDistY < sideDistZ)
		{
			sideDistY += deltaDY;
			mapY += stepY;
			side = 1;
		} 
		
		else 
		{
			sideDistZ += deltaDZ;
			mapZ += stepZ;
			side = 2;
		}

		float block = GetVoxel(ivec3(mapX, mapY, mapZ));

		if (block != 0) 
		{
			hit = 1;
			blockType = block;

			if (side == 0) 
			{
				T = (mapX - origin.x + (1 - stepX) / 2) / direction.x + t1;
				normal = vec3(1, 0, 0) * -stepX;
			}

			else if (side == 1) 
			{
				T = (mapY - origin.y + (1 - stepY) / 2) / direction.y + t1;
				normal = vec3(0, 1, 0) * -stepY;
			}

			else
			{
				T = (mapZ - origin.z + (1 - stepZ) / 2) / direction.z + t1;
				normal = vec3(0, 0, 1) * -stepZ;
			}

			break;
		}
	}

	return T;
}

float raySphereIntersect(vec3 r0, vec3 rd, vec3 s0, float sr) 
{
    float a = dot(rd, rd);
    vec3 s0_r0 = r0 - s0;
    float b = 2.0 * dot(rd, s0_r0);
    float c = dot(s0_r0, s0_r0) - (sr * sr);

    if ( b * b - 4.0 * a * c < 0.0) 
    {
        return -1.0;
    }

    return (-b - sqrt((b * b) - 4.0 * a * c)) / (2.0 * a);
}

void main()
{
	// checker boarding test
	//if (int(gl_FragCoord.x + gl_FragCoord.y) % 2 == (u_CurrentFrame % 2))
	//{
	//	o_Color = vec3(0.0f);
	//	return;
	//}

    Ray r;
    r.Origin = v_RayOrigin;
    r.Direction = normalize(v_RayDirection);

	vec3 normal;
	float id;
	float t = voxel_traversal(r.Origin, r.Direction, normal, id);
	bool intersect = t > 0.0f;
    vec3 world_position = intersect ? r.Origin + (r.Direction * t) : vec3(-1.0f);
	vec2 UV;
	int normal_idx = 0;

	CalculateUV(world_position, normal, UV, normal_idx);

	if (intersect)
	{
		o_Position.xyz = world_position;
		o_Normal = normal;
	} 

	else
	{
		o_Position.xyz = vec3(-1.0f);
		o_Normal = vec3(-1.0f);
	}

	o_Position.w = t;

	int reference_id = clamp(int(floor(id * 255.0f)), 0, 127);
	vec3 texture_ids = BLOCK_TEXTURE_DATA[reference_id];

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

	o_Data = vec4(texture_ids.x, texture_ids.y, texture_ids.z, 1.0f);
}