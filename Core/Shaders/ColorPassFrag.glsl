#version 330 core

#define PCF_COUNT 6
#define SMOOTH_SHADOW_SAMPLING

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;
in vec3 v_RayDirection;

uniform sampler2D u_DiffuseTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InitialTracePositionTexture;
uniform sampler2D u_DataTexture;
uniform sampler2D u_ShadowTexture;
uniform sampler2DArray u_BlockAlbedoTextures;
uniform sampler2DArray u_BlockNormalTextures;
uniform sampler2DArray u_BlockPBRTextures;
uniform sampler2DArray u_BlueNoiseTextures;
uniform samplerCube u_Skybox;

uniform vec3 SunDirection;
uniform vec3 MoonDirection;

uniform mat4 u_ShadowView;
uniform mat4 u_ShadowProjection;

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

vec4 ClampedTexture(sampler2D tex, vec2 txc)
{
    return texture(tex, clamp(txc, 0.0f, 1.0f));
}

vec3 NeighbourhoodClamping(vec3 tempColor, sampler2D tex, vec2 txc) 
{
	vec2 neighbourhoodOffsets[8] = vec2[8]
	(
		vec2(-1.0, -1.0),
		vec2( 0.0, -1.0),
		vec2( 1.0, -1.0),
		vec2(-1.0,  0.0),
		vec2( 1.0,  0.0),
		vec2(-1.0,  1.0),
		vec2( 0.0,  1.0),
		vec2( 1.0,  1.0)
	);

	vec3 minclr = vec3(0.0f), maxclr = vec3(1.0);
    vec2 View = 1.0f / textureSize(tex, 0);

	for(int i = 0; i < 8; i++) 
	{
		vec2 offset = neighbourhoodOffsets[i] * View;
		vec3 clr = texture(tex, txc + offset, 0.0).rgb;
		minclr = min(minclr, clr);
		maxclr = max(maxclr, clr);
	}

	return clamp(tempColor, minclr, maxclr);
}

vec3 BilateralUpsample(sampler2D tex, vec2 txc, vec3 base_normal, float base_depth)
{
    const vec2 Kernel[4] = vec2[](
        vec2(0.0f, 1.0f),
        vec2(1.0f, 0.0f),
        vec2(-1.0f, 0.0f),
        vec2(0.0, -1.0f)
    );

    vec2 texel_size = 1.0f / textureSize(tex, 0);

    vec3 color = vec3(0.0f, 0.0f, 0.0f);
    float weight_sum;

    for (int i = 0; i < 4; i++) 
    {
        vec3 sampled_normal = texture(u_NormalTexture, txc + Kernel[i] * texel_size).xyz;
        float nweight = pow(abs(dot(sampled_normal, base_normal)), 32);

        float sampled_depth = texture(u_InitialTracePositionTexture, txc + Kernel[i] * texel_size).z; 
        float dweight = 1.0f / (abs(base_depth - sampled_depth) + 0.001f);

        float computed_weight = nweight * dweight;
        color.rgb += texture(tex, txc + Kernel[i] * texel_size).rgb * computed_weight;
        weight_sum += computed_weight;
    }

    color /= max(weight_sum, 0.2f);
    color = clamp(color, texture(tex, txc).rgb * 0.4f, vec3(1.0f));
    return color;
}

const vec2 PoissonDisk[32] = vec2[]
(
    vec2(-0.613392, 0.617481),  vec2(0.751946, 0.453352),
    vec2(0.170019, -0.040254),  vec2(0.078707, -0.715323),
    vec2(-0.299417, 0.791925),  vec2(-0.075838, -0.529344),
    vec2(0.645680, 0.493210),   vec2(0.724479, -0.580798),
    vec2(-0.651784, 0.717887),  vec2(0.222999, -0.215125),
    vec2(0.421003, 0.027070),   vec2(-0.467574, -0.405438),
    vec2(-0.817194, -0.271096), vec2(-0.248268, -0.814753),
    vec2(-0.705374, -0.668203), vec2(0.354411, -0.887570),
    vec2(0.977050, -0.108615),  vec2(0.175817, 0.382366),
    vec2(0.063326, 0.142369),   vec2(0.487472, -0.063082),
    vec2(0.203528, 0.214331),   vec2(-0.084078, 0.898312),
    vec2(-0.667531, 0.326090),  vec2(0.488876, -0.783441),
    vec2(-0.098422, -0.295755), vec2(0.470016, 0.217933),
    vec2(-0.885922, 0.215369),  vec2(-0.696890, -0.549791),
    vec2(0.566637, 0.605213),   vec2(-0.149693, 0.605762),
    vec2(0.039766, -0.396100),  vec2(0.034211, 0.979980)
);

vec2 ReprojectShadow(in vec3 world_pos)
{
	vec3 WorldPos = world_pos;

	vec4 ProjectedPosition = u_ShadowProjection * u_ShadowView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;

	return ProjectedPosition.xy;
}

float ComputeShadow(in vec2 txc)
{
    float shadow;

#ifdef SMOOTH_SHADOW_SAMPLING
    vec2 TexSize = textureSize(u_ShadowTexture, 0);
    vec2 TexelSize = 1.0 / TexSize; 

	for(int x = 0; x <= PCF_COUNT; x++)
	{
        float noise = texture(u_BlueNoiseTextures, vec3(gl_FragCoord.xy / textureSize(u_BlueNoiseTextures, 0).xy, 0.0f)).r;
        float theta = noise * 6.28318530718;
        float cosTheta = cos(theta);
        float sinTheta = sin(theta);
        mat2 dither = mat2(vec2(cosTheta, -sinTheta), vec2(sinTheta, cosTheta));

		vec2 jitter_value;
        jitter_value = PoissonDisk[x] * dither;

        float pcf = texture(u_ShadowTexture, txc + jitter_value * TexelSize).r; 
		shadow += pcf;        
	}

	shadow /= float(PCF_COUNT);
#else 
    shadow = textureBicubic(u_ShadowTexture, txc).r;
#endif

    return shadow;
}

bool IsInScreenSpaceBounds(in vec2 tx)
{
    if (tx.x > 0.0f && tx.y > 0.0f && tx.x < 1.0f && tx.y < 1.0f)
    {
        return true;
    }

    return false;
}

void main()
{
    vec4 WorldPosition = texture(u_InitialTracePositionTexture, v_TexCoords);
    vec3 SampledNormals = texture(u_NormalTexture, v_TexCoords).rgb;

    o_Color = vec3(1.0f);
    o_Color = texture(u_Skybox, normalize(v_RayDirection)).rgb * 0.5f;

    if (WorldPosition.w > 0.0f)
    {
        vec2 uv;
        vec2 ReprojectedShadowPos = ReprojectShadow(WorldPosition.xyz);
        float RayTracedShadow = 1.0f;
        
        if (IsInScreenSpaceBounds(ReprojectedShadowPos))
        {
            RayTracedShadow = 1.0f - ComputeShadow(ReprojectedShadowPos);
        }

        CalculateUV(WorldPosition.xyz, SampledNormals, uv); uv.y = 1.0f - uv.y;

        vec2 TexSize = textureSize(u_InitialTracePositionTexture, 0);
        float PixelDepth1 = texture(u_InitialTracePositionTexture, clamp(v_TexCoords + vec2(0.0f, 1.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;
        float PixelDepth2 = texture(u_InitialTracePositionTexture, clamp(v_TexCoords + vec2(0.0f, -1.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;
        float PixelDepth3 = texture(u_InitialTracePositionTexture, clamp(v_TexCoords + vec2(1.0f, 0.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;
        float PixelDepth4 = texture(u_InitialTracePositionTexture, clamp(v_TexCoords + vec2(-1.0f, 0.0f) * (1.0f / TexSize), 0.001f, 0.999f)).w;

        if (PixelDepth1 > 0.0f && PixelDepth2 > 0.0f && PixelDepth3 > 0.0f && PixelDepth4 > 0.0f)
        {
            vec4 data = texture(u_DataTexture, v_TexCoords);
            vec3 AlbedoColor = texture(u_BlockAlbedoTextures, vec3(uv, data.x)).rgb;
            vec3 NormalMap = texture(u_BlockNormalTextures, vec3(uv, data.y)).rgb;
            vec3 PBRMap = texture(u_BlockPBRTextures, vec3(uv, data.z)).rgb;

            //vec3 Diffuse = BilateralUpsample(u_DiffuseTexture, v_TexCoords, SampledNormals, WorldPosition.z).rgb;
            vec3 Diffuse = textureBicubic(u_DiffuseTexture, v_TexCoords).rgb;

            vec3 LightAmbience = (vec3(120.0f, 172.0f, 255.0f) / 255.0f) * 1.01f;
            vec3 Ambient = (AlbedoColor * LightAmbience) * 0.005;

            o_Color = Ambient + ((AlbedoColor * 1.7f) * (Diffuse * 1.2f)) * RayTracedShadow;
            return;
        }
    }

    else 
    {
        vec3 ray_dir = normalize(v_RayDirection); bool intersect_body;

        if(dot(ray_dir, normalize(SunDirection)) > 0.9997f)
        {
            o_Color = vec3(4.0f) * 3.0f; intersect_body = true;
        }

        if(dot(ray_dir, normalize(MoonDirection)) > 0.99986f)
        {
            o_Color = vec3(0.6f, 0.6f, 0.9f) * 50.0f; intersect_body = true;
        }
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
