#version 330 core

#define CHUNK_SIZE_X 64
#define CHUNK_SIZE_Y 64
#define CHUNK_SIZE_Z 64

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

float GetVoxel(vec3 loc)
{
    if (IsInVoxelizationVolume(loc))
    {
         return texture(u_VoxelDataTexture, loc).r;
    }
    
    return 0.0f;
}

float GetVoxel(ivec3 loc)
{
    if (IsInVoxelizationVolume(loc))
    {
         return texture(u_VoxelDataTexture, vec3(loc)).r;
    }
    
    return 0.0f;
}

bool VoxelExists(in vec3 loc)
{
    if (GetVoxel(loc) > 0.0f) 
    {
        return true;
    }

    return false;
}

// Calculates uv from world position
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv)
{
    world_pos.x = abs(world_pos.x);
    world_pos.y = abs(world_pos.y);
    world_pos.z = abs(world_pos.z);

    const vec3 NORMAL_TOP = vec3(0.0f, 1.0f, 0.0f);
    const vec3 NORMAL_BOTTOM = vec3(0.0f, -1.0f, 0.0f);
    const vec3 NORMAL_FRONT = vec3(0.0f, 0.0f, 1.0f);
    const vec3 NORMAL_BACK = vec3(0.0f, 0.0f, -1.0f);
    const vec3 NORMAL_LEFT = vec3(-1.0f, 0.0f, 0.0f);
    const vec3 NORMAL_RIGHT = vec3(1.0f, 0.0f, 0.0f);

    if (normal == NORMAL_TOP)
    {
        uv = vec2(fract(world_pos.xz));
    }

    else if (normal == NORMAL_BOTTOM)
    {
        uv = vec2(fract(world_pos.zx));
    }

    else if (normal == NORMAL_RIGHT)
    {
        uv = vec2(fract(world_pos.xy));
    }

    else if (normal == NORMAL_LEFT)
    {
        uv = vec2(fract(world_pos.yx));
    }
    
    else if (normal == NORMAL_FRONT)
    {
        uv = vec2(fract(world_pos.zy));
    }

     else if (normal == NORMAL_BOTTOM)
    {
        uv = vec2(fract(world_pos.yz));
    }
}

// By Snurf (Thanks snurf!)
bool RaytraceVoxel(Ray r, out float voxel, out vec3 hitPos, out ivec3 hitIndex, out vec3 hitNormal, out vec2 uv, const int dist) 
{
    vec3 origin = r.Origin;
    vec3 direction = r.Direction;

    ivec3 mapPos = ivec3(floor(origin));
    vec3 deltaDist = abs(vec3(length(direction)) / direction);
    ivec3 rayStep = ivec3(sign(direction));
    vec3 sideDist = (sign(direction) * (vec3(mapPos) - origin) + (sign(direction) * 0.5) + 0.5) * deltaDist; 

    bvec3 mask = bvec3(false);
    bool hit = false;
    float id = -1;
    
    for (int i = 0; i < dist && !hit; i++) 
    {
        voxel = GetVoxel(mapPos);
        id = voxel;

        if (id > 0)
        {
            hitPos = direction / sum(vec3(mask) * direction) * sum(vec3(mask) * (mapPos + vec3(lessThan(direction, vec3(0))) - origin)) + origin;
            hit = true;
        }

        mask = lessThanEqual(sideDist.xyz, min(sideDist.yzx, sideDist.zxy));
        sideDist += vec3(mask) * deltaDist;
        mapPos += ivec3(vec3(mask)) * rayStep;

        hitNormal = rayStep * vec3(mask);
        CalculateUV((hitPos), vec3(mask), uv);
    }
    
    hitIndex = mapPos;

    return hit;
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
    //if (int(gl_FragCoord.x + gl_FragCoord.y) % 2 == 0)
    //{
    //    discard;
    //}
    //
    //Ray r;
    //r.Origin = v_RayOrigin;
    //r.Direction = normalize(v_RayDirection);
    //
    //float voxel;
    //vec3 hitpos;
    //ivec3 hitidx;
    //vec3 normal;
    //vec2 uv = vec2(1.0f);
    //
    //bool hit_voxel = RaytraceVoxel(r, voxel, hitpos, hitidx, normal, uv, 256);
    //
    //if (hit_voxel)
    //{
    //    o_Color = vec3(uv, 0.0f);
    //    return;
    //}
    //
	//o_Color = GetSkyColorAt(r.Direction);

    o_Color = vec3(texture(u_VoxelDataTexture, vec3(5, 5, 5)).r);
}