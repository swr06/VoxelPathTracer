#version 330 core

layout (location = 0) out float o_OutputShadow;

uniform sampler2D u_Texture;
uniform sampler2D u_BlockIDs;
uniform sampler2D u_Depth;

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

float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))

#define PI 3.141592653

vec4 SampleShadow(sampler2D tex, vec2 TexCoord) {
    if (TexCoord != clamp(TexCoord, 0.001f, 0.999f)) {
        return vec4(0.0f);
    }

    return texture(tex, clamp(TexCoord,0.000001f,0.999999f));
}

const vec2 poisson60[60] = vec2[60]( vec2( 0.0f, 0.0f ), vec2( -0.6882537922853785f, -0.7111317767221583f ), vec2( -0.6274096173225843f, 0.7656851631168115f ), vec2( 0.438935753135918f, -0.8979943375165644f ), vec2( 0.4502162526854448f, 0.8548678722681659f ), vec2( 0.9845638474773585f, -0.07095114459273748f ), vec2( -0.8934233017286533f, 0.03908749241352094f ), vec2( -0.1041298091721986f, -0.5798950770680016f ), vec2( -0.06803356559479536f, 0.6102496213286652f ), vec2( 0.509174241198879f, 0.2894054214168671f ), vec2( 0.38689048423442296f, -0.3738296618653449f ), vec2( -0.46071005827985306f, 0.11656661447622403f ), vec2( 0.7630966823560734f, -0.5656398395450662f ), vec2( -0.9172362926277846f, -0.37535593154110897f ), vec2( -0.48677392389692337f, -0.37385894915219947f ), vec2( -0.3031348988374572f, -0.9466992193949866f ), vec2( 0.7859202342028092f, 0.5606120923356787f ), vec2( -0.867211464270054f, 0.44597828133370276f ), vec2( 0.0323423512036535f, 0.9794273613211496f ), vec2( 0.70357087824837f, -0.2671229361213631f ), vec2( 0.027544262276607227f, -0.976726632912478f ), vec2( -0.4560087702465609f, 0.4663848193859052f ), vec2( 0.1239400070732809f, 0.35343694593150926f ), vec2( 0.9023373945449029f, 0.255114354359261f ), vec2( 0.3440127031811375f, -0.04121145996349828f ), vec2( -0.29915370950834713f, 0.8654864939577471f ), vec2( -0.2411135572484363f, -0.1770342883360836f ), vec2( 0.09525157553205667f, -0.32545303446259854f ), vec2( 0.53766788718549f, 0.595960436402185f ), vec2( 0.2617483050036872f, -0.6649446275742013f ), vec2( -0.1932151444452412f, 0.22207797827158834f ), vec2( -0.37949874264978245f, -0.664944787001864f ), vec2( 0.6958027736103634f, 0.02793970483741243f ), vec2( -0.6926665307044337f, -0.172050962314387f ), vec2( -0.7156188590221898f, 0.24024198757139947f ), vec2( 0.1900140953184063f, 0.6931433340653482f ), vec2( 0.9292678293470402f, -0.3582768846319978f ), vec2( 0.5584488633584952f, -0.7028393286643093f ), vec2( -0.13384215441365507f, -0.3636780197733723f ), vec2( 0.3331912420320544f, 0.4632002058767273f ), vec2( 0.23250403745290502f, 0.15219064833479282f ), vec2( -0.4211413452338232f, 0.6841988195095107f ), vec2( 0.014067119190510312f, -0.766973667257993f ), vec2( -0.667440732080303f, 0.03612363793081464f ), vec2( -0.2448993011248424f, 0.47582438419349743f ), vec2( -0.9821464192949045f, -0.16396565723304637f ), vec2( -0.47840420460294264f, -0.14610676387278682f ), vec2( -0.5830666284680038f, -0.5411700088274732f ), vec2( -0.657946007912339f, 0.5235313850402264f ), vec2( -0.5052662734490019f, -0.8395535562823743f ), vec2( 0.6994981412038591f, 0.35281977338813797f ), vec2( -0.8170741335631682f, -0.5590451991736097f ), vec2( 0.574922793347461f, -0.12579738935151363f ), vec2( 0.56367944237224f, -0.47988494783886176f ), vec2( 0.23946039004927042f, 0.9434345037259955f ), vec2( 0.2425494748436006f, -0.8680924838483116f ), vec2( 0.5156169195725088f, 0.08185666945403365f ), vec2( 0.004384418000440115f, 0.7927069742243914f ), vec2( 0.22790860430470603f, -0.48480025691005396f ), vec2( -0.20664196341432908f, 0.022917405996018454f ) );

void main()
{
    float id = texelFetch(u_BlockIDs, ivec2(v_TexCoords * textureSize(u_BlockIDs, 0).xy), 0).r;
	int iid = clamp(int(floor(id * 255.0f)), 0, 127);
    int SSSMaskFetch = BlockSSSData[iid];

    if (SSSMaskFetch > 0 && v_TexCoords == clamp(v_TexCoords, 0.04f, 0.96f)) {

        float TotalShadow = 0.0f;
        vec2 TexelSize = 1.0f / textureSize(u_Texture, 0).xy;
        float LinearDepth = texture(u_Depth, v_TexCoords).x;
        float Hash = bayer32(gl_FragCoord.xy);
        float Theta = Hash * 2.0f * PI;
        float CosTheta = cos(Theta);
        float SinTheta = sin(Theta);
        mat2 RotationMatrix = mat2(vec2(CosTheta, -SinTheta), vec2(SinTheta, CosTheta));
        const float Diagonal = sqrt(2.0f);
        const int Samples = 32;
		int SuccessfulSamples = 0;

        for (int x = 0; x < Samples ; x++) {
            vec2 SampleCoord = v_TexCoords + (RotationMatrix * poisson60[x]) * TexelSize * 16.0f;
            
            if (SampleCoord == clamp(SampleCoord, 0.0f, 1.0)) 
			{
                //int BlockIDAt = clamp(int(floor(( texelFetch(u_BlockIDs, ivec2(SampleCoord * textureSize(u_BlockIDs, 0).xy), 0).r) * 255.0f)), 0, 127);
                //if (BlockIDAt == iid || BlockIDAt == 0)
                
                float DepthAt = texture(u_Depth, SampleCoord).x;

                bool Skybias = DepthAt < 0.0f; //&& distance(SampleCoord, v_TexCoords) < 1. / 8.;

                if (abs(DepthAt - LinearDepth) <= Diagonal || Skybias) 
                {
                    TotalShadow += SampleShadow(u_Texture, SampleCoord).x;
					SuccessfulSamples++;
                }
            }
        }

        TotalShadow /= float(Samples);
        o_OutputShadow = TotalShadow;
		
		if (SuccessfulSamples <= 1) {
			float BaseShadow = texture(u_Texture,v_TexCoords).x;
			o_OutputShadow = BaseShadow;
		}
    }

    else {
        float BaseShadow = texture(u_Texture,v_TexCoords).x;
        o_OutputShadow = BaseShadow;
    }
}