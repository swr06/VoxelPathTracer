// Ground fog ->

#version 330 core


layout (location = 0) out vec3 o_Volumetrics;


uniform sampler2D u_Mask;
uniform samplerCube u_Skymap;
uniform sampler2D u_GBuffer;

uniform vec3 u_ViewerPosition;
uniform vec2 u_Dimensions;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_Time;
uniform float u_Strength;



// Noise FBMs 
// Perlin noise
float mod289(float x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 mod289(vec4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 perm(vec4 x){return mod289(((x * 34.0) + 1.0) * x);}

float noise(vec3 p)
{
    vec3 a = floor(p);
    vec3 d = p - a;
    d = d * d * (3.0 - 2.0 * d);
    vec4 b = a.xxyy + vec4(0.0, 1.0, 0.0, 1.0);
    vec4 k1 = perm(b.xyxy);
    vec4 k2 = perm(k1.xyxy + b.zzww);
    vec4 c = k2 + a.zzzz;
    vec4 k3 = perm(c);
    vec4 k4 = perm(c + 1.0);
    vec4 o1 = fract(k3 * (1.0 / 41.0));
    vec4 o2 = fract(k4 * (1.0 / 41.0));
    vec4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    vec2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);
    return o4.y * d.y + o4.x * (1.0 - d.y);
}


// Simplex noise 
vec3 random3(vec3 c) 
{
	float j = 4096.0*sin(dot(c,vec3(17.0, 59.4, 15.0)));
	vec3 r;
	r.z = fract(512.0*j);
	j *= .125;
	r.x = fract(512.0*j);
	j *= .125;
	r.y = fract(512.0*j);
	return r-0.5;
}

float simplex3d(vec3 p) {
    const float F3 =  0.3333333f;
    const float G3 =  0.1666667f;
    vec3 s = floor(p + dot(p, vec3(F3)));
    vec3 x = p - s + dot(s, vec3(G3));
    vec3 e = step(vec3(0.0), x - x.yzx);
    vec3 i1 = e*(1.0 - e.zxy);
    vec3 i2 = 1.0 - e.zxy*(1.0 - e);
    vec3 x1 = x - i1 + G3;
    vec3 x2 = x - i2 + 2.0*G3;
    vec3 x3 = x - 1.0 + 3.0*G3;
    vec4 w, d;
    w.x = dot(x, x);
    w.y = dot(x1, x1);
    w.z = dot(x2, x2);
    w.w = dot(x3, x3);
    w = max(0.6 - w, 0.0);
    d.x = dot(random3(s), x);
    d.y = dot(random3(s + i1), x1);
    d.z = dot(random3(s + i2), x2);
    d.w = dot(random3(s + 1.0), x3);
    w *= w;
    w *= w;
    d *= w;
    return dot(d, vec4(52.0));
}

const mat3 RotationMatrix_1 = mat3(-0.37f, 0.36f, 0.85f,-0.14f,-0.93f, 0.34f,0.92f, 0.01f,0.40f);
const mat3 RotationMatrix_2 = mat3(-0.55f,-0.39f, 0.74f, 0.33f,-0.91f,-0.24f,0.77f, 0.12f,0.63f);
const mat3 RotationMatrix_3 = mat3(-0.71f, 0.52f,-0.47f,-0.08f,-0.72f,-0.68f,-0.7f,-0.45f,0.56f);

// Rotate each sample octave to reduce artifacts 
// Works super well usually 
float SimplexFractalNoise(vec3 m) {
    return  4.0f *(  0.5333333f * simplex3d(m * RotationMatrix_1)
			+ 0.2666667f * simplex3d(2.0f * m * RotationMatrix_2)
			+ 0.1333333f * simplex3d(4.0f * m * RotationMatrix_3)
			+ 0.0666667f * simplex3d(8.0f * m));
}


float ScatterIntegral(float OD, float Coefficent)
{
    float InverseLog2 = 1.0f / log(2.0);
    float T = -Coefficent * InverseLog2; 
    float B = -1.0f / Coefficent;
    float C =  1.0f / Coefficent;
    //return exp(T * OD) * B + C;
    return exp2(T * OD) * B + C;
}


// Returns the optical depth for a point using a perlin/simplex fractal FBM
float OpticalDepth(vec3 p)
{
    if (true) {

        //float OD = pow(dot(noise(p - Time),vec3(0.625,0.25,0.625),2.2)*4.0 ;
        //OD -= pow(dot(noise(p - Time*2.0f),vec3(0.625,0.25,0.625),8.2)*4.0 ;
        //OD += pow(dot(noise(p + Time*4.0f),vec3(0.625,0.25,0.625),4.2)*4.0 ;
        //float Integral = ScatterIntegral(OD,1.1f);

        float Time = mod(u_Time * (0.9f / 1.0f), 100.0f);
        float OD = noise(p - Time) * 0.5f; 
        p *= 1.2f;
        OD += noise(p * 2.0f + Time) * 0.25f;
        OD += noise(p * 4.0f - Time) * 0.125f;
        float OpticalDepthSquared = OD * OD;
	    return (OpticalDepthSquared * 4.0f + 0.25f) * (0.9f / 1.0f);
    }

    else {
        float Time = u_Time * 0.9f;
        Time = mod(Time, 100.0f);
        p*=3.0f;
        float OD = SimplexFractalNoise(vec3(p.x - mod(Time, 100.0f), p.y + Time * 0.2f, p.z + Time)) * 1.0f;
	    return (OD*OD * 4.0f + 0.25f) * (0.9f / 1.0f) * (0.9f / 1.0f);
    }
}

void main() 
{


}