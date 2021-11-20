// shadertoy my beloved

#define HASHSCALE1 .1031
#define HASHSCALE3 vec3(.1031, .1030, .0973)
#define HASHSCALE4 vec4(1031, .1030, .0973, .1099)

float hash11(float p) {
    vec3 p3  = fract(vec3(p) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

float hash12(vec2 p) {
    vec3 p3  = fract(vec3(p.xyx) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

float hash13(vec3 p3) {
    p3  = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

vec2 hash22(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.xx + p3.yz) * p3.zy);

}

vec2 hash23(vec3 p3) {
    p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.xx + p3.yz) * p3.zy);
}

vec3 hash32(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yxz + 19.19);
    return fract((p3.xxy + p3.yzz) * p3.zyx);
}

vec3 hash33(vec3 p3) {
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yxz + 19.19);
    return fract((p3.xxy + p3.yxx) * p3.zyx);
}

vec4 hash42(vec2 p) {
	vec4 p4 = fract(vec4(p.xyxy) * HASHSCALE4);
    p4 += dot(p4, p4.wzxy + 19.19);
    return fract((p4.xxyz + p4.yzzw) * p4.zywx);

}

vec4 hash43(vec3 p) {
	vec4 p4 = fract(vec4(p.xyzx)  * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

uint initRand(uint val0, uint val1, uint backoff)
{
    backoff = 16;
	uint v0 = val0, v1 = val1, s0 = 0;

	for (uint n = 0; n < backoff; n++)
	{
		s0 += 0x9e3779b9;
		v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s0) ^ ((v1 >> 5) + 0xc8013ea4);
		v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s0) ^ ((v0 >> 5) + 0x7e95761e);
	}
	return v0;
}

// Takes our seed, updates it, and returns a pseudorandom float in [0..1]
float nextRand(inout uint s)
{
	s = (1664525u * s + 1013904223u);
	return float(s & 0x00FFFFFF) / float(0x01000000);
}

vec4 noiseSmooth(vec2 coord) {
    coord = coord * noiseTextureResolution + 0.5;

	vec2 whole = floor(coord);
	vec2 part  = cubicSmooth(fract(coord));

	coord = (whole + part - 0.5) * noiseResInverse;

	return texture2D(noisetex, coord);
}

float Get3DNoise(vec3 pos) {
	float p = floor(pos.z);
	float f = pos.z - p;
	
	const float zStretch = 17.0 * noiseResInverse;
	
	vec2 coord = pos.xy * noiseResInverse + (p * zStretch);
	
	vec2 noise = noiseSmooth(coord).xy;
	
	return mix(noise.x, noise.y, f);
}

float worley12(vec2 n) {
    float dis = 1.0;
    for(int x = -1; x <= 1; x++) {
        for(int y = -1; y <= 1; y++) {
            vec2  p = floor(n) + vec2(x, y);
            float d = length(hash22(p) + vec2(x, y) - fract(n));
            if (dis>d) dis = d;
        }
    }
    return clamp01(1.0 - dis);
}

float worley13(vec3 n) 
{
    float dis = 1.0;
    for(int x = -1; x <= 1; x++) {
        for(int y = -1; y <= 1; y++) {
            for (int z = -1; z <= 1; z++) {
                vec3  p = floor(n) + vec3(x, y, z);
                float d = length(hash33(p) + vec3(x, y, z) - fract(n));
                if (dis>d) dis = d;
            }
        }
    }
    return clamp01(1.0 - dis);
}

vec3 Permutation(vec3 x) 
{
  return mod((34.0 * x + 1.0) * x, 289.0);
}

vec2 inoise(vec3 P, float jitter)
{			
	vec3 Pi = mod(floor(P), 289.0);
 	vec3 Pf = fract(P);
	vec3 oi = vec3(-1.0, 0.0, 1.0);
	vec3 of = vec3(-0.5, 0.5, 1.5);
	vec3 px = Permutation(Pi.x + oi);
	vec3 py = Permutation(Pi.y + oi);

	vec3 p, ox, oy, oz, dx, dy, dz;
	vec2 F = vec2(1e6);

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			p  = Permutation(px[i] + py[j] + Pi.z + oi); // pij1, pij2, pij3

			ox = fract(p*K) - Ko;
			oy = mod(floor(p*K),7.0)*K - Ko;
			
			p  = Permutation(p);
			
			oz = fract(p*K) - Ko;
		
			dx = Pf.x - of[i] + jitter*ox;
			dy = Pf.y - of[j] + jitter*oy;
			dz = Pf.z - of + jitter*oz;
			
			vec3 d = dx * dx + dy * dy + dz * dz; // dij1, dij2 and dij3, squared
			
			//Find lowest and second lowest distances
			for(int n = 0; n < 3; n++)
			{
				if(d[n] < F[0])
				{
					F[1] = F[0];
					F[0] = d[n];
				}
				else if(d[n] < F[1])
				{
					F[1] = d[n];
				}
			}
		}
	}
	
	return F;
}