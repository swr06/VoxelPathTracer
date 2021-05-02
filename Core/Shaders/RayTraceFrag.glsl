#version 330 core

#define CHUNK_SIZE_X 1024
#define CHUNK_SIZE_Y 128
#define CHUNK_SIZE_Z 1024

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

uniform sampler3D u_VoxelDataTexture;

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};


vec3 MapSize = vec3(CHUNK_SIZE_X, CHUNK_SIZE_Y, CHUNK_SIZE_Z);


vec3 GetSkyColorAt(vec3 rd) 
{
    vec3 unit_direction = normalize(rd);

    float t = 0.5f * (unit_direction.y + 1.0);
    return (1.0 - t) * vec3(1.0, 1.0, 1.0) +  t * vec3(0.5, 0.7, 1.0);
}

float sum(vec3 v)
{
    return v.x + v.y + v.z;
}

bool IsInVoxelizationVolume(in vec3 pos)
{
    if (pos.x < 0.0f || pos.y < 0.0f || pos.z < 0.0f || 
        pos.x > float(CHUNK_SIZE_X) || pos.y > float(CHUNK_SIZE_Y) || pos.z > float(CHUNK_SIZE_Z))
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
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv)
{
    const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
    const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
    const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
    const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
    const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
    const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

    if (normal == NORMAL_TOP)
    {
        uv = vec2(mod(world_pos.xz, 1.0f));
    }

    else if (normal == NORMAL_BOTTOM)
    {
        uv = vec2(mod(world_pos.xz, 1.0f));
    }

    else if (normal == NORMAL_RIGHT)
    {
        uv = vec2(mod(world_pos.zy, 1.0f));
    }

    else if (normal == NORMAL_LEFT)
    {
        uv = vec2(mod(world_pos.zy, 1.0f));
    }
    
    else if (normal == NORMAL_FRONT)
    {
        uv = vec2(mod(world_pos.xy, 1.0f));
    }

     else if (normal == NORMAL_BACK)
    {
        uv = vec2(mod(world_pos.xy, 1.0f));
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

	for (int i = 0; i < 400; i++) 
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
    Ray r;
    r.Origin = v_RayOrigin;
    r.Direction = normalize(v_RayDirection);
    
	vec3 normal;
	float id;

	float t = voxel_traversal(r.Origin, r.Direction, normal, id);
    vec3 world_position = r.Origin + (r.Direction * t);
	bool intersect = t > 0.0f;
	vec2 UV;

	CalculateUV(world_position, normal, UV);

    if (intersect)
    {
        o_Color = vec3(UV, 0.0f);
        return;
    }
    
	o_Color = GetSkyColorAt(r.Direction);
}