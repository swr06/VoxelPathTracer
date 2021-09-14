float ValueNoise(vec3 pos) 
{
    vec3 p  = floor(pos); 
    vec3 b  = fract(pos);

    vec2 uv = (p.xy + (vec2(-97.0) * p.z)) + b.xy;
    vec2 rg = texture(noisetex, (uv)/256.0).xy;

    return cubic_smooth(mix(rg.x, rg.y, b.z));
}

vec3 Get2DNoiseFromTex(vec2 pos) {
    return texture(noisetex, pos).xyz;
}

//x-large scale worley, y-small scale worley, z-perlin-worley
vec3 Get3DNoiseFromTex(vec3 pos) 
{   
	return texture(depthtex2, fract(pos)).xyz;
}

//x-large scale worley, y-small scale worley, z-perlin-worley
vec3 Get3DNoiseFromTex(vec2 coord) 
{   
	vec3 pos    = vec3(coord.x, 0.0, coord.y);
	return noise_3d(pos);
}