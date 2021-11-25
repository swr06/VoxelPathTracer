#version 330 core
#define R_INNER 0.985f

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_ViewPosition;
in vec3 v_RayDirection;

uniform float u_Time;
uniform vec3 u_SunDirection;

float MiePhase(float g, float c, float cc) 
{
	float gg = g * g;
	float a = (1.0f - gg) * (1.0f + cc);
	float b = 1.0f + gg - 2.0 * g * c;
	b *= sqrt(b);
	b *= 2.0f + gg;	
	return 1.5f * a / b;
}

float ReyleighPhase(float cc) 
{
	return 0.75f * (1.0f + cc);
}

vec2 RSI(vec3 p, vec3 dir, float r) 
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

float SampleDensity(vec3 p)
{
	const float R = 1.0;
	const float SCALE_H = 4.0 / (R - R_INNER);
	const float SCALE_L = 1.0 / (R - R_INNER);
	return exp(-(length(p) - R_INNER) * SCALE_H) * 2.0f;
}

float SampleOpticalDepth(vec3 p, vec3 q) 
{
	const int Steps = 3;

	vec3 Increment = (q - p) / float(Steps);
	Increment /= 4.0f;
	vec3 v = p + Increment * 0.5f;
	float Accumulated = 0.0f;

	for (int i = 0; i < Steps; i++) 
	{
		Accumulated += SampleDensity(v);
		v += Increment;
	}

	float StepSize = length(Increment);
	Accumulated *= StepSize * (1.0 / (1.0f - R_INNER));
	return Accumulated;
}

vec3 PrimaryScatter(vec3 Origin, vec3 Direction, vec2 e, vec3 l, const float MieAmount, const float RayleighAmount) 
{
	const float PI = 3.14159265359f;
	const float R = 1.0f;
	const float SCALE_L = 1.0f / (R - R_INNER);
	const vec3 RayleighScatterCoefficient = vec3(0.2f, 0.45f, 1.0f);	
	const float G_M = -0.75f;

	const float Steps = 12;
	const float K_R = 0.186f * RayleighAmount;
	const float K_M = 0.035f * MieAmount;
	
	// Fake boost :
	float boosty = clamp((l.y + 0.1f), 0.0f, 1.0f) * 0.95f + 0.05f;
	boosty = 1.0f / sin(boosty);
	float Length = (e.y * (1.0 + boosty * 0.0)) / float(Steps);
	
	vec3 Increment = Direction * Length;
	Increment *= 2.0f;

	vec3 RayPosition = Origin;
	vec3 v = RayPosition + Direction * (Length * (0.5f + boosty * 0.0f));

	vec3 Accumulated = vec3(0.0f);

	for (int i = 0; i < Steps; i++) 
	{
		vec2 f = RSI(v, l, R);
		vec3 u = v + l * f.y;
		float n = (SampleOpticalDepth(RayPosition, v) + SampleOpticalDepth(v, u)) * (PI * 4.0f);
		Accumulated += SampleDensity(v) * exp(-n * (K_R * RayleighScatterCoefficient + K_M));
		v += Increment;
	}

	Accumulated *= Length * SCALE_L;

	float c = dot(Direction, -l);
	float cc = c * c;
	vec2 Phase = vec2(ReyleighPhase(cc),MiePhase(G_M, c, cc));
	return Accumulated * (K_R * RayleighScatterCoefficient * Phase.x + K_M * Phase.y) * 13.0f;
}

vec3 SecondScatter(vec3 o, vec3 dir, vec2 e, vec3 l) 
{
	const float Steps = 8;
	const float PI = 3.14159265359;
	const float R = 1.0;
	const float SCALE_L = 1.0 / (R - R_INNER);
	const float K_R = 0.166;
	const float K_M = 0.00;
	const float E = 14.3;
	const vec3 C_R = vec3(0.2, 0.6, 1.0);	
	const float G_M = -0.65;

	float len = (e.y) / float(Steps);
	vec3 step = dir * len;
	step *= 2.0;
	vec3 p = o;

	// Fake boost :
	float boosty = clamp((l.y + 0.1), 0.0f, 1.0f) * 0.95 + 0.05;
	boosty = 1.0 / sin(boosty);


	vec3 v = p + dir * (len * (0.5 + boosty * 0.0));
	vec3 sum = vec3(0.0);

	for (int i = 0; i < Steps; i++) 
	{
		vec2 f = RSI(v, l, R);
		vec3 u = v + l * f.y;
		float n = (SampleOpticalDepth(p, v) + SampleOpticalDepth(v, u)) * (PI * 4.0);
		sum += SampleDensity(v) * exp(-n * (K_R * C_R + K_M));
		v += step;
	}


	sum *= len * SCALE_L;
	float c = dot(dir, -l);
	float cc = c * c;

	vec2 Phase = vec2(ReyleighPhase(cc), MiePhase(G_M, c, cc)); 
	return sum * (K_R * C_R * Phase.x + K_M * Phase.y) * E;
}

vec3 AtmosphericScattering(vec3 rayDir, vec3 lightVector, const float mieAmount)
{
	// Constants :
	const float PI = 3.14159265359;
	const float DEG_TO_RAD = PI / 180.0;
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

	vec3 Eye = vec3(0.0, mix(R_INNER, 1.0, 0.05), 0.0);
	vec3 OriginalRayDirection = rayDir;

	// Fake boost :
	float boosty = clamp((lightVector.y), 0.0f, 1.0f) * 0.90 + 0.10;
	boosty = 1.0 / sin(boosty);

	if (rayDir.y < 0.0)
	{
		rayDir.y = abs(rayDir.y);
		rayDir.y *= rayDir.y;
	}

	vec3 Up = vec3(0.0f, 1.0f, 0.0f);
	vec2 IntersectedEye = RSI(Eye, rayDir, R);
	vec2 IntersectedUp = RSI(Eye, Up, R);

	// March!
	vec3 MarchedAtmosphere = PrimaryScatter(Eye, rayDir, IntersectedEye, lightVector, mieAmount, 1.0);
	vec3 Secondary = SecondScatter(Eye, Up, IntersectedUp, lightVector);
	
	vec3 AmbientTerm = vec3(0.3f, 0.5f, 1.0f);

	MarchedAtmosphere += dot(Secondary, vec3(0.66)) * AmbientTerm;
	MarchedAtmosphere *= vec3(0.8, 0.89, 1.0);
	MarchedAtmosphere = pow(MarchedAtmosphere, vec3(1.2));


	return MarchedAtmosphere;
}

void main()
{    
	vec3 col = AtmosphericScattering(normalize(v_RayDirection), normalize(u_SunDirection), 0.12f);
	vec3 col2 = AtmosphericScattering(normalize(v_RayDirection), normalize(vec3(-u_SunDirection.x, -u_SunDirection.y, u_SunDirection.z)), 0.5f);
    float SunVisibility = clamp(dot(u_SunDirection, vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0; SunVisibility = 1.0f  - SunVisibility;
	o_Color = vec4(mix(col*0.6f*vec3(1.), col2*((vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.05f), SunVisibility), 1.0f);
}