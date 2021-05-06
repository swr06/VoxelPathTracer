#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;

uniform sampler2D u_DiffuseTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InitialTracePositionTexture;
uniform sampler2D u_DataTexture;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform samplerCube u_Skybox;

vec4 textureBicubic(sampler2D sampler, vec2 texCoords);
void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv);

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

//Due to low sample count we "tonemap" the inputs to preserve colors and smoother edges
vec3 WeightedSample(sampler2D colorTex, vec2 texcoord)
{
	vec3 wsample = texture(colorTex,texcoord).rgb * 1.0f;
	return wsample / (1.0f + GetLuminance(wsample));
}

vec3 smoothfilter(in sampler2D tex, in vec2 uv)
{
	vec2 textureResolution = textureSize(tex, 0);
	uv = uv*textureResolution + 0.5;
	vec2 iuv = floor( uv );
	vec2 fuv = fract( uv );
	uv = iuv + fuv*fuv*fuv*(fuv*(fuv*6.0-15.0)+10.0);
	uv = (uv - 0.5)/textureResolution;
	return WeightedSample( tex, uv);
}

vec3 sharpen(in sampler2D tex, in vec2 coords) 
{
	vec2 renderSize = textureSize(tex, 0);
	float dx = 1.0 / renderSize.x;
	float dy = 1.0 / renderSize.y;
	vec3 sum = vec3(0.0);
	sum += -1. * smoothfilter(tex, coords + vec2( -1.0 * dx , 0.0 * dy));
	sum += -1. * smoothfilter(tex, coords + vec2( 0.0 * dx , -1.0 * dy));
	sum += 5. * smoothfilter(tex, coords + vec2( 0.0 * dx , 0.0 * dy));
	sum += -1. * smoothfilter(tex, coords + vec2( 0.0 * dx , 1.0 * dy));
	sum += -1. * smoothfilter(tex, coords + vec2( 1.0 * dx , 0.0 * dy));
	return sum;
}


void main()
{
    vec4 WorldPosition = texture(u_InitialTracePositionTexture, v_TexCoords);
    vec3 SampledNormals = texture(u_NormalTexture, v_TexCoords).rgb;

    o_Color = vec3(1.0f);

    if (WorldPosition.w > 0.0f)
    {
        float id = texture(u_DataTexture, v_TexCoords).r;
        vec3 Diffuse = WeightedSample(u_DiffuseTexture, v_TexCoords).rgb;

        vec2 uv;
        CalculateUV(WorldPosition.xyz, SampledNormals, uv);

        //o_Color = vec3(uv, 1.0f);
        o_Color = vec3(texture(u_BlockAlbedoTextures, vec3(uv, id))) * Diffuse;
    }

    else 
    {
        o_Color = texture(u_Skybox, normalize(v_RayDirection)).rgb;
    }
}

vec4 cubic(float v){
    vec4 n = vec4(1.0, 2.0, 3.0, 4.0) - v;
    vec4 s = n * n * n;
    float x = s.x;
    float y = s.y - 4.0 * s.x;
    float z = s.z - 4.0 * s.y + 6.0 * s.x;
    float w = 6.0 - x - y - z;
    return vec4(x, y, z, w) * (1.0/6.0);
}

vec4 textureBicubic(sampler2D sampler, vec2 texCoords)
{

   vec2 texSize = textureSize(sampler, 0);
   vec2 invTexSize = 1.0 / texSize;

   texCoords = texCoords * texSize - 0.5;


    vec2 fxy = fract(texCoords);
    texCoords -= fxy;

    vec4 xcubic = cubic(fxy.x);
    vec4 ycubic = cubic(fxy.y);

    vec4 c = texCoords.xxyy + vec2 (-0.5, +1.5).xyxy;

    vec4 s = vec4(xcubic.xz + xcubic.yw, ycubic.xz + ycubic.yw);
    vec4 offset = c + vec4 (xcubic.yw, ycubic.yw) / s;

    offset *= invTexSize.xxyy;

    vec4 sample0 = texture(sampler, offset.xz);
    vec4 sample1 = texture(sampler, offset.yz);
    vec4 sample2 = texture(sampler, offset.xw);
    vec4 sample3 = texture(sampler, offset.yw);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix(
       mix(sample3, sample2, sx), mix(sample1, sample0, sx)
    , sy);
}

void CalculateUV(vec3 world_pos, in vec3 normal, out vec2 uv)
{
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
        uv = vec2(fract(world_pos.xz));
    }

    else if (normal == NORMAL_RIGHT)
    {
        uv = vec2(fract(world_pos.zy));
    }

    else if (normal == NORMAL_LEFT)
    {
        uv = vec2(fract(world_pos.zy));
    }
    
    else if (normal == NORMAL_FRONT)
    {
        uv = vec2(fract(world_pos.xy));
    }

     else if (normal == NORMAL_BACK)
    {
        uv = vec2(fract(world_pos.xy));
    }
}
