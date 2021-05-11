#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384

#define ALBEDO_TEX_LOD 3 // 512, 256, 128

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InitialTraceNormalTexture;
uniform sampler2D u_PBRTexture;

uniform sampler3D u_VoxelData;
uniform samplerCube u_Skymap;

uniform vec3 BLOCK_TEXTURE_DATA[128];
uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockAlbedoTextures;

uniform vec3 u_ViewerPosition;
		
// Function prototypes
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv);
void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv);
float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType, in int mdist) ;
float GetVoxel(ivec3 loc);
bool IsInVoxelizationVolume(in vec3 pos);

void main()
{
	// Start ray at sampled position, use normalized normal (already in tangent space) as direction, trace and get the albedo color at.
	
	vec3 SampledWorldPosition = texture(u_PositionTexture, v_TexCoords).rgb;
	vec3 NormalMap = texture(u_NormalTexture, v_TexCoords).rgb;
	vec3 InitialTraceNormal = texture(u_InitialTraceNormalTexture, v_TexCoords).rgb;

	vec3 NormalMapped = normalize(InitialTraceNormal + (NormalMap * 0.25f));
	float MetalnessAt = texture(u_PBRTexture, v_TexCoords).g;

	if (MetalnessAt > 0.1f)
	{
		vec3 I = normalize(SampledWorldPosition - u_ViewerPosition);
		vec3 R = normalize(reflect(I, NormalMapped));
		vec3 Normal;
		float Blocktype;

		float T = voxel_traversal(SampledWorldPosition, R, Normal, Blocktype, 100);
		vec3 HitPosition = SampledWorldPosition + (R * T);

		if (T > 0.0f)
		{
			int reference_id = clamp(int(floor(Blocktype * 255.0f)), 0, 127);
			vec3 texture_ids = BLOCK_TEXTURE_DATA[reference_id];

			vec2 UV; 
			CalculateUV(HitPosition, Normal, UV);

			o_Color = textureLod(u_BlockAlbedoTextures, vec3(UV,texture_ids.x), ALBEDO_TEX_LOD).rgb;
		}

		else 
		{
			o_Color = texture(u_Skymap, R).rgb;
		}
	}

	else 
	{
		o_Color = vec3(-1.0f);
	}

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
         return texelFetch(u_VoxelData, loc, 0).r;
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

float ProjectToCube(vec3 ro, vec3 rd) 
{	
	const vec3 MapSize = vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);
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

float voxel_traversal(vec3 orig, vec3 direction, inout vec3 normal, inout float blockType, in int mdist) 
{
	const vec3 MapSize = vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);

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

	for (int i = 0; i < mdist; i++) 
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

void CalculateVectors(vec3 world_pos, in vec3 normal, out vec3 tangent, out vec3 bitangent, out vec2 uv)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );

	const vec3 Tangents[6] = vec3[]( vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f),
					 vec3(0.0f, 0.0f, -1.0f), vec3(0.0f, 0.0f, -1.0f)
				   );

	const vec3 BiTangents[6] = vec3[]( vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f),
				     vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 1.0f),
					 vec3(0.0f, -1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f)
	);

	if (normal == Normals[0])
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[0];
		bitangent = BiTangents[0];
    }

    else if (normal == Normals[1])
    {
        uv = vec2(fract(world_pos.xy));
		tangent = Tangents[1];
		bitangent = BiTangents[1];
    }

    else if (normal == Normals[2])
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[2];
		bitangent = BiTangents[2];
    }

    else if (normal == Normals[3])
    {
        uv = vec2(fract(world_pos.xz));
		tangent = Tangents[3];
		bitangent = BiTangents[3];
    }
	
    else if (normal == Normals[4])
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[4];
		bitangent = BiTangents[4];
    }
    

    else if (normal == Normals[5])
    {
        uv = vec2(fract(world_pos.zy));
		tangent = Tangents[5];
		bitangent = BiTangents[5];
    }
}

void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv)
{
	// Hard coded normals, tangents and bitangents

    const vec3 Normals[6] = vec3[]( vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, -1.0f),
					vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, -1.0f, 0.0f), 
					vec3(-1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f)
			      );

	if (normal == Normals[0])
    {
        uv = vec2(fract(world_pos.xy));
    }

    else if (normal == Normals[1])
    {
        uv = vec2(fract(world_pos.xy));
    }

    else if (normal == Normals[2])
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (normal == Normals[3])
    {
        uv = vec2(fract(world_pos.xz));
    }
	
    else if (normal == Normals[4])
    {
        uv = vec2(fract(world_pos.zy));
    }
    

    else if (normal == Normals[5])
    {
        uv = vec2(fract(world_pos.zy));
    }
}