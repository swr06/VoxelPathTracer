#version 330 core
#define PI 3.1415926536

/*
https://www.shadertoy.com/view/XtBXDz
*/

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_ViewPosition;
in vec3 v_RayDirection;

uniform float u_Time;
uniform vec3 u_SunDirection;

uniform int u_NumSamples;
uniform int u_NumLightSamples; 

// Structures
struct ray_t
{
    vec3 origin;
    vec3 direction;
};

struct sphere_t
{
    vec3 origin;
    float radius;
};
    
const vec3 betaR = vec3(5.5e-6, 13.0e-6, 22.4e-6); // Rayleigh 
const vec3 betaM = vec3(21e-6); // Mie
const float hR = 7994.0; // Rayleigh
const float hM = 1200.0; // Mie
float g_Time = 0.0f;

float rayleigh_phase_func(float mu)
{
    return 3.0f * (1.0f + mu * mu) / (16.0f * PI);
}

//A beloved classic
const float g = 0.76;
float henyey_greenstein_phase_func(float mu)
{
	return
						(1. - g*g)
	/ //---------------------------------------------
		((4. * PI) * pow(1. + g*g - 2.*g*mu, 1.5));
}

bool intersects_sphere(in ray_t ray, in sphere_t sphere, inout float t0, inout float t1)
{
    vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	float d2 = dot(rc, rc) - tca * tca;
	if (d2 > radius2) return false;
	float thc = sqrt(radius2 - d2);
	t0 = tca - thc;
	t1 = tca + thc;
	return true;
}


const float earth_radius = 6360e3; 
const float atmosphere_radius = 6420e3;

const vec3 sun_power = vec3(16.0) * vec3(0.8f, 0.85f, 1.0f);
const vec3 moon_power = vec3(4.0) * ((vec3(96.0f, 192.0f, 255.0f) / 255.0f) * 0.6f);

const sphere_t atmosphere = sphere_t(vec3(0.),atmosphere_radius);
const sphere_t earth = sphere_t(vec3(0.),earth_radius);

bool get_sun_light(in ray_t ray, inout float opticalDepthR, inout float opticalDepthM)
{
    float t0,t1;
    intersects_sphere(ray, atmosphere,t0,t1);
    
    float march_pos = 0.0f ;
    float march_length = t1 / float(u_NumLightSamples);
    
    for(int i = 0; i < u_NumLightSamples; i++)
    {
        vec3 s = ray.origin + ray.direction * (march_pos + 0.5f * march_length);
        float height = length(s)-earth_radius;
        opticalDepthR += exp(-height/hR) * march_length;
        opticalDepthM += exp(-height/hM) * march_length;
        
        march_pos += march_length;
    }

    return true;
}

vec3 get_incident_light(in ray_t ray)
{
    vec3 sun_dir = u_SunDirection ; 
    vec3 moon_dir = vec3(-sun_dir.x, -sun_dir.y, sun_dir.z); 
    
    float t0,t1;

    if(!intersects_sphere(ray, atmosphere, t0, t1))
    {
        return vec3(1.0f);
    }
    
    float march_length = t1 / float(u_NumSamples);
    
    float mu = dot(sun_dir, ray.direction);
    float muMoon = dot(moon_dir, ray.direction);
    
    float phaseR = rayleigh_phase_func(mu) * 1.66f;
    float phaseM = henyey_greenstein_phase_func(mu) * 1.05f;
    float phaseMoonR = rayleigh_phase_func(muMoon);
    float phaseMoonM = henyey_greenstein_phase_func(muMoon) * 1.0;
    float opticalDepthR = 0.0f;
    float opticalDepthM = 0.0f;
    float opticalDepthMoonR = 0.0f;
    float opticalDepthMoonM = 0.0f;
    
    vec3 sumR = vec3(0.0f);
    vec3 sumM = vec3(0.0f);
    vec3 sumMoonM = vec3(0.0f);
    vec3 sumMoonR = vec3(0.0f);
    
    float march_pos = 0.;
    
    for(int i = 0 ; i < u_NumSamples; i++)
    {
        vec3 s = ray.origin + ray.direction * (march_pos + 0.5f * march_length); //sample middle of step
        float height = length(s) - earth_radius;
        
        float hr = exp(-height / hR) * march_length;
        float hm = exp(-height / hM) * march_length;
        opticalDepthR += hr;
        opticalDepthM += hm;
        
        ray_t light_ray = ray_t(s, sun_dir);
        ray_t moon_ray = ray_t(s,moon_dir);
        
        float opticalDepthLightR = 0.;
        float opticalDepthLightM = 0.;
        float opticalDepthLightMoonR = 0.;
        float opticalDepthLightMoonM = 0.;
        
        bool overground = get_sun_light(light_ray, opticalDepthLightR,  opticalDepthLightM);
        bool overgroundMoon = get_sun_light(moon_ray, opticalDepthLightMoonR, opticalDepthLightMoonM);
        
        if(overground || overgroundMoon)
        {
            vec3 tau = betaR * (opticalDepthR + opticalDepthLightR) 
                + betaM * 1.1 * (opticalDepthM + opticalDepthLightM);
            
            vec3 tauMoon = betaR * (opticalDepthR + opticalDepthLightMoonR) 
                + betaM * 1.1 * (opticalDepthM + opticalDepthLightMoonM);

            vec3 attenuation = exp(-tau);
            vec3 attenuationMoon = exp(-tauMoon);
            
            sumR += hr * attenuation;
            sumM += hm * attenuation;
            sumMoonR += hr * attenuationMoon;
            sumMoonM += hm * attenuationMoon;
        }

        march_pos += march_length;
    }
    
    vec3 col = sun_power * (sumR * phaseR * betaR +sumM * phaseM * betaM) 
        + moon_power * (sumMoonR * phaseMoonR * betaR + sumMoonM * phaseMoonM * betaM); 
    
    return col;
}

void main()
{    
    g_Time = u_Time * 0.1f;
    vec2 uv = v_TexCoords.xy;
    vec3 ray_dir = v_RayDirection;
    ray_dir.y = clamp(ray_dir.y, 0.02f, 1.0f);

    vec3 cameraCenter = vec3(0.0f, earth_radius + 1.0f, 0.0f);
    ray_t primary_ray = ray_t(cameraCenter, normalize(ray_dir)); 
    vec3 col = get_incident_light(primary_ray);
    col = 1.0 - exp(-col * 0.96f);
    o_Color = vec4(col, 1.0f);
}