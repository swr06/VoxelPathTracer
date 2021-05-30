#version 330 core

/*
Traversal Paper used : https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
*/

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

#define MULTIPLE_TEXTURING_GRASS
#define ALPHA_TESTING

layout (location = 0) out vec4 o_Position;
layout (location = 1) out vec3 o_Normal;
layout (location = 2) out vec4 o_Data;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform int u_CurrentFrame;

uniform sampler3D u_VoxelDataTexture;
uniform sampler2DArray u_AlbedoTextures;

uniform vec2 u_Dimensions;
uniform vec4 BLOCK_TEXTURE_DATA[128];
uniform float BLOCK_EMISSIVE_TEXTURE_DATA[128];

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
const vec3 NORMALS[6] = vec3[]
(
					vec3(1.0, 0.0, 0.0),
					vec3(-1.0, 0.0, 0.0),
					vec3(0.0, 1.0, 0.0),
					vec3(0.0, -1.0, 0.0),
					vec3(0.0, 0.0, 1.0),
					vec3(0.0, 0.0, -1.0)
);

float sum(vec3 v)
{
    return v.x + v.y + v.z;
}

bool InInVoxelVolume(in vec3 pos)
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
    if (InInVoxelVolume(loc))
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

//float ProjectToCube(vec3 ro, vec3 rd) 
//{	
//	float tx1 = (0 - ro.x) / rd.x;
//	float tx2 = (MapSize.x - ro.x) / rd.x;
//
//	float ty1 = (0 - ro.y) / rd.y;
//	float ty2 = (MapSize.y - ro.y) / rd.y;
//
//	float tz1 = (0 - ro.z) / rd.z;
//	float tz2 = (MapSize.z - ro.z) / rd.z;
//
//	float tx = max(min(tx1, tx2), 0);
//	float ty = max(min(ty1, ty2), 0);
//	float tz = max(min(tz1, tz2), 0);
//
//	float t = max(tx, max(ty, tz));
//	
//	return t;
//}

bool voxel_traversal(Ray r, inout float block, out vec3 normal, out vec3 world_pos)
{
	world_pos = r.Origin;
	vec2 txc;

	vec3 Temp;
	vec3 VoxelCoord; 
	vec3 FractPosition;

	Temp.x = r.Direction.x > 0.0 ? 1.0 : 0.0;
	Temp.y = r.Direction.y > 0.0 ? 1.0 : 0.0;
	Temp.z = r.Direction.z > 0.0 ? 1.0 : 0.0;

	vec3 plane = floor(world_pos + Temp);

	for (int x = 0; x < 225; x++)
	{
		if (!InInVoxelVolume(world_pos))
		{
			break;
		}

		vec3 Next = (plane - world_pos) / r.Direction;
		int side = 0;

		if(x > 0) 
		{
			if (Next.x < min(Next.y, Next.z)) 
			{
				world_pos += r.Direction * Next.x;
				world_pos.x = plane.x;
				plane.x += sign(r.Direction.x);
				side = 0;
			}

			else if (Next.y < Next.z) 
			{
				world_pos += r.Direction * Next.y;
				world_pos.y = plane.y;
				plane.y += sign(r.Direction.y);
				side = 1;
			}

			else 
			{
				world_pos += r.Direction * Next.z;
				world_pos.z = plane.z;
				plane.z += sign(r.Direction.z);
				side = 2;
			}
		}

		VoxelCoord = (plane - Temp);
		FractPosition = fract(world_pos);

		switch (side)
		{
			case 0:
			{
				txc = FractPosition.zy;
				break;
			}

			case 1:
			{
				txc = FractPosition.xz;
				break;
			}

			default:
			{
				txc = FractPosition.xy;
				break;
			}
		}

		txc.y = 1.0f - txc.y;

		int Side = ((side + 1) * 2) - 1;

		if (side == 0) 
		{
			if (world_pos.x - VoxelCoord.x > 0.5)
			{
				Side = 0;
			}
		}

		else if (side == 1)
		{
			if (world_pos.y - VoxelCoord.y > 0.5)
			{
				Side = 2;
			}
		}

		else 
		{
			if (world_pos.z - VoxelCoord.z > 0.5)
			{
				Side = 4;
			}
		}

		normal = NORMALS[Side];
		block = GetVoxel(ivec3(VoxelCoord.xyz));

		#ifdef ALPHA_TESTING

		if (block > 0)
		{
			int reference_id = clamp(int(floor(block * 255.0f)), 0, 127);
			bool transparent = BLOCK_TEXTURE_DATA[reference_id].a > 0.5f;

			if (transparent)
			{
				int temp_idx; 

				if (texture(u_AlbedoTextures, vec3(vec2(txc.x, txc.y), BLOCK_TEXTURE_DATA[reference_id].x)).a < 0.1f)
				{
					continue;
				}
			}

			return true; 
		}
		#endif
	}

	return false;
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
    Ray r;
    r.Origin = v_RayOrigin;
    r.Direction = normalize(v_RayDirection);

	int normal_idx = 0;

	vec3 normal;
	float id;
	int face;
	vec2 UV;

	vec3 world_position;
	bool intersect = voxel_traversal(r, id, normal, world_position);
	float t = intersect ? distance(r.Origin, world_position) : -1.0f;

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

	int reference_id;
	vec4 texture_ids;
	bool transparent;

	if (intersect)
	{
		reference_id = clamp(int(floor(id * 255.0f)), 0, 127);
		texture_ids.xyz = BLOCK_TEXTURE_DATA[reference_id].rgb;
		texture_ids.w = BLOCK_EMISSIVE_TEXTURE_DATA[reference_id];
	}

	else 
	{
		texture_ids = vec4(-1.0f);
	}

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