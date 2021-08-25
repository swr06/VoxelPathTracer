#version 330 core

#define WORLD_SIZE_X 384
#define WORLD_SIZE_Y 128
#define WORLD_SIZE_Z 384
#define PI 3.14159265359

#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))


layout (location = 0) out vec3 o_Volumetrics;

in vec2 v_TexCoords;

uniform sampler3D u_ParticipatingMedia;

uniform vec3 u_ViewerPosition;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

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
		vec3 float_loc = vec3(loc) / vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z); // normalize
        return texture(u_ParticipatingMedia, float_loc).r;
    }
    
    return 0.0f;
}

float saturate(float x) {
	return clamp(x, 1e-5f, 1.0f);
}

void main() 
{
	// Ray properties
	vec3 rO = u_ViewerPosition;

	float Dither = bayer64(gl_FragCoord.xy);
	vec3 rD = GetRayDirectionAt(v_TexCoords);
	rD = normalize(rD);


	vec3 direction = rD;

	float TotalDensity = 0.0f;

	// Utility for DDA 
	vec3 world_pos = rO;
	vec3 Temp;
	vec3 VoxelCoord; 
	vec3 FractPosition;
	Temp.x = rD.x > 0.0 ? 1.0 : 0.0;
	Temp.y = rD.y > 0.0 ? 1.0 : 0.0;
	Temp.z = rD.z > 0.0 ? 1.0 : 0.0;
	vec3 plane = floor(world_pos + Temp);


	for (int x = 0; x < 30; x++)
	{
		if (!IsInVolume(world_pos))
		{
			break;
		}

		vec3 Next = (plane - world_pos) / direction;
		int side = 0;

		if (Next.x < min(Next.y, Next.z)) 
		{
			world_pos += direction * Next.x;
			world_pos.x = plane.x;
			plane.x += sign(direction.x);
			side = 0;
		}

		else if (Next.y < Next.z) 
		{
			world_pos += direction * Next.y;
			world_pos.y = plane.y;
			plane.y += sign(direction.y);
			side = 1;
		}

		else 
		{
			world_pos += direction * Next.z;
			world_pos.z = plane.z;
			plane.z += sign(direction.z);
			side = 2;
		}

		VoxelCoord = (plane - Temp);

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

		TotalDensity += saturate(GetVoxel(ivec3(VoxelCoord.xyz)) * 3.0f );
	}

	o_Volumetrics = vec3(TotalDensity);
}