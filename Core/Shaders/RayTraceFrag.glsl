#version 330 core

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;
in vec3 v_RayOrigin;

struct Ray
{
	vec3 Origin;
	vec3 Direction;
};

float raySphereIntersect(vec3 r0, vec3 rd, vec3 s0, float sr) 
{
    float a = dot(rd, rd);
    vec3 s0_r0 = r0 - s0;
    float b = 2.0 * dot(rd, s0_r0);
    float c = dot(s0_r0, s0_r0) - (sr * sr);

    if (b*b - 4.0*a*c < 0.0) 
    {
        return -1.0;
    }

    return (-b - sqrt((b*b) - 4.0*a*c))/(2.0*a);
}

void main()
{
    Ray r;
    r.Origin = v_RayOrigin;
    r.Direction = normalize(v_RayDirection);

    if (raySphereIntersect(r.Origin, r.Direction, vec3(0.0f, 0.0f, -1.0f), 0.5f) > 0)
    {
    	o_Color = vec4(1.0f, 0.0f, 0.0f, 1.0f);
        return;
    }

	o_Color = vec4(1.0f, 1.0f, 1.0f, 1.0f);
}