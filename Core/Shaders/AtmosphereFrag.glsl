#version 330 core
#define R_INNER 0.985f

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_ViewPosition;
in vec3 v_RayDirection;

uniform float u_Time;
uniform vec3 u_SunDirection;

float phase_mie(float g, float c, float cc) 
{
	float gg = g * g;
	float a = (1.0 - gg) * (1.0 + cc);
	float b = 1.0 + gg - 2.0 * g * c;
	b *= sqrt( b );
	b *= 2.0 + gg;	
	return 1.5 * a / b;
}

float phase_reyleigh(float cc) 
{
	return 0.75 * (1.0 + cc);
}

vec2 RaySphereIntersection( vec3 p, vec3 dir, float r ) 
{
	float b = dot(p, dir);
	float c = dot(p, p) - r * r;
	
	float d = b * b - c;
	if (d < 0.0) 
	{
		return vec2(10000.0, -10000.0);
	}

	d = sqrt(d);
	
	return vec2( -b - d, -b + d );
}

float density(vec3 p)
{
	const float R = 1.0;
	const float SCALE_H = 4.0 / (R - R_INNER);
	const float SCALE_L = 1.0 / (R - R_INNER);

	return exp(-(length(p) - R_INNER) * SCALE_H) * 2.0;
}

float optic(vec3 p, vec3 q) 
{
	const int numOutscatter = 1;

	const float R = 1.0;
	const float SCALE_L = 1.0 / (R - R_INNER);

	vec3 step = (q - p) / float(numOutscatter);
	step *= 0.3;

	vec3 v = p + step * 0.5;
	float sum = 0.0;

	for (int i = 0; i < numOutscatter; i++) 
	{
		sum += density(v);
		v += step;
	}

	sum *= length(step) * SCALE_L;
	return sum;
}

vec3 in_scatter(vec3 o, vec3 dir, vec2 e, vec3 l, const float mieAmount, const float rayleighAmount) 
{
	const float numInscatter = 4;
	
	const float PI = 3.14159265359;

	const float R = 1.0;
	const float SCALE_L = 1.0 / (R - R_INNER);

	const float K_R = 0.186 * rayleighAmount;
	const float K_M = 0.035 * mieAmount;
	const float E = 14.3;
	const vec3 C_R = vec3(0.2, 0.45, 1.0);	//Rayleigh scattering coefficients
	const float G_M = -0.75;

	float boosty = clamp((l.y + 0.1), 0.0f, 1.0f) * 0.95 + 0.05;
	boosty = 1.0 / sin(boosty);

	float len = (e.y * (1.0 + boosty * 0.0)) / float(numInscatter);
	vec3 step = dir * len;
	step *= 2.0;
	vec3 p = o;

	vec3 v = p + dir * (len * (0.5 + boosty * 0.0));
	vec3 sum = vec3(0.0);

	for ( int i = 0; i < numInscatter; i++ ) 
	{
		vec2 f = RaySphereIntersection( v, l, R );
		vec3 u = v + l * f.y;
		
		float n = (optic(p, v) + optic(v, u)) * (PI * 4.0);
		sum += density(v) * exp(-n * (K_R * C_R + K_M));
		v += step;
	}

	sum *= len * SCALE_L;
	float c  = dot(dir, -l);
	float cc = c * c;
	return sum * (K_R * C_R * phase_reyleigh(cc) + K_M * phase_mie(G_M, c, cc)) * E;
}

vec3 in_scatter2(vec3 o, vec3 dir, vec2 e, vec3 l) 
{
	const float numInscatter = 2;
	
	const float PI = 3.14159265359;

	const float R = 1.0;
	const float SCALE_L = 1.0 / (R - R_INNER);

	const float K_R = 0.166;
	const float K_M = 0.00;
	const float E = 14.3;
	const vec3 C_R = vec3(0.2, 0.6, 1.0);	//Rayleigh scattering coefficients
	const float G_M = -0.65;

	float len = (e.y) / float(numInscatter);
	vec3 step = dir * len;
	step *= 2.0;
	vec3 p = o;

	//float boosty = 1.0 - abs(l.y);
	float boosty = clamp((l.y + 0.1), 0.0f, 1.0f) * 0.95 + 0.05;
	boosty = 1.0 / sin(boosty);
	vec3 v = p + dir * (len * (0.5 + boosty * 0.0));
	vec3 sum = vec3(0.0);

	for (int i = 0; i < numInscatter; i++) 
	{
		vec2 f = RaySphereIntersection(v, l, R);
		vec3 u = v + l * f.y;
		float n = (optic(p, v) + optic(v, u)) * (PI * 4.0);
		sum += density(v) * exp(-n * (K_R * C_R + K_M));
		v += step;
	}

	sum *= len * SCALE_L;
	float c  = dot(dir, -l);
	float cc = c * c;
	return sum * (K_R * C_R * phase_reyleigh(cc) + K_M * phase_mie(G_M, c, cc)) * E;
}

vec3 AtmosphericScattering(vec3 rayDir, vec3 lightVector, const float mieAmount)
{
	const float PI = 3.14159265359;
	const float DEG_TO_RAD = PI / 180.0;

	//Scatter constants
	const float K_R = 0.166;
	const float K_M = 0.0025;
	const float E = 14.3;
	const vec3 C_R = vec3(0.3, 0.7, 1.0);	//Rayleigh scattering coefficients
	const float G_M = -0.85;

	const float R = 1.0;
	const float SCALE_H = 4.0 / (R - R_INNER);
	const float SCALE_L = 1.0 / (R - R_INNER);

	const int NUM_OUT_SCATTER = 10;
	const float FNUM_OUT_SCATTER = 10.0;

	const int NUM_IN_SCATTER = 10;
	const float FNUM_IN_SCATTER = 10.0;

	vec3 eye = vec3(0.0, mix(R_INNER, 1.0, 0.05), 0.0);

	vec3 originalRayDir = rayDir;

	if (rayDir.y < 0.0)
	{
		rayDir.y = abs(rayDir.y);
		rayDir.y *= rayDir.y;
	}

	vec3 up = vec3(0.0, 1.0, 0.0);
	vec2 e = RaySphereIntersection(eye, rayDir, R);
	vec2 eup = RaySphereIntersection(eye, up, R);
	vec3 atmosphere = in_scatter(eye, rayDir, e, lightVector, mieAmount, 1.0);
	vec3 secondary = in_scatter2(eye, up, eup, lightVector);
	vec3 ambient = vec3(0.3, 0.5, 1.0);
	vec3 ground = vec3(0.1, 0.1, 0.1) * 0.05;
	float boosty = clamp((lightVector.y), 0.0f, 1.0f) * 0.90 + 0.10;
	boosty = 1.0 / sin(boosty);

	atmosphere += dot(secondary, vec3(0.66)) * ambient;
	atmosphere *= vec3(0.8, 0.89, 1.0);
	atmosphere = pow(atmosphere, vec3(1.2));
	return atmosphere;
}

void main()
{    
	vec3 col = AtmosphericScattering(normalize(v_RayDirection), normalize(u_SunDirection), 0.12f);
	vec3 col2 = AtmosphericScattering(normalize(v_RayDirection), normalize(vec3(-u_SunDirection.x, -u_SunDirection.y, u_SunDirection.z)), 0.5f);
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	o_Color = vec4(mix(col*0.625f, col2*((vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.06f), SunVisibility), 1.0f);
}