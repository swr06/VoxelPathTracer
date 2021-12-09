#version 330 core

layout (location = 0) out float o_OutputShadow;

uniform sampler2D u_Texture;
uniform sampler2D u_BlockIDs;

in vec2 v_TexCoords;

layout (std430, binding = 0) buffer SSBO_BlockData
{
    int BlockAlbedoData[128];
    int BlockNormalData[128];
    int BlockPBRData[128];
    int BlockEmissiveData[128];
	int BlockTransparentData[128];
	int BlockSSSData[128];
};

vec4 SampleShadow(vec2 TexCoord) {
    if (TexCoord != clamp(TexCoord, 0.001f, 0.999f)) {
        return vec4(0.0f);
    }

    return texture(u_Texture, clamp(TexCoord,0.000001f,0.999999f));
}

float GaussianVertical(vec2 uv) 
{
    vec2 Resolution = textureSize(u_Texture, 0);
    vec2 Direction = vec2(0.0f, 1.0f);
    float Shadow = 0.0f;
    vec2 Offset = vec2(1.333333) * Direction;
    Shadow += SampleShadow(uv).x * 0.294117;
    Shadow += SampleShadow(uv + (Offset / Resolution)).x * 0.352941;
    Shadow += SampleShadow(uv - (Offset / Resolution)).x * 0.352941;
    return Shadow; 
}
  
float GaussianShadow(vec2 uv) 
{
    vec2 Resolution = textureSize(u_Texture, 0);
    vec2 Direction = vec2(1.0f, 0.0f);
    float Shadow = 0.0f;
    vec2 Offset = vec2(1.333333) * Direction;
    Shadow += GaussianVertical(uv) * 0.294117f;
    Shadow += GaussianVertical(uv + (Offset / Resolution)) * 0.352941f;
    Shadow += GaussianVertical(uv - (Offset / Resolution)) * 0.352941f;
    return Shadow; 
}

void main()
{
    float id = texelFetch(u_BlockIDs, ivec2(v_TexCoords * textureSize(u_BlockIDs, 0).xy), 0).r;
	int iid = clamp(int(floor(id * 255.0f)), 0, 127);
    int SSSFetch = BlockSSSData[iid];

    if (SSSFetch > 0) {
        o_OutputShadow = GaussianShadow(v_TexCoords);
    }

    else {
        float BaseShadow = texture(u_Texture,v_TexCoords).x;
        o_OutputShadow = BaseShadow;
    }
}