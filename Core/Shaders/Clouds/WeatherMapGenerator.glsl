#version 450 core

// a 1024x1024 texture

layout (location = 0) out vec4 o_Noise;

in vec2 v_TexCoords;

vec2 seed = vec2(2.0f, 3.9684f);

float saturate(float x) {
    return clamp(x,0.0f,1.0f);
}

float random2D( in vec2 st ) 
{
    return fract( sin(dot( st.xy, vec2(12.9898,78.233 ) + seed.xy ) ) * 43758.5453123);
}

uniform float perlinAmplitude = 0.5;
uniform float perlinFrequency = 0.8;
uniform float perlinScale = 100.0;
uniform int perlinOctaves = 4;

float noiseInterpolation(in vec2 i_coord, in float i_size)
{
	vec2 grid = i_coord * i_size;
    
    vec2 randomInput = floor( grid );
    vec2 weights     = fract( grid );
    
    
    float p0 = random2D( randomInput );
    float p1 = random2D( randomInput + vec2( 1.0, 0.0  ) );
    float p2 = random2D( randomInput + vec2( 0.0, 1.0 ) );
    float p3 = random2D( randomInput + vec2( 1.0, 1.0 ) );
    
    weights = smoothstep( vec2( 0.0, 0.0 ), vec2( 1.0, 1.0 ), weights ); 
    
    return p0 +
           ( p1 - p0 ) * ( weights.x ) +
           ( p2 - p0 ) * ( weights.y ) * ( 1.0 - weights.x ) +
           ( p3 - p1 ) * ( weights.y * weights.x );    
}

float perlinNoise(vec2 uv, float sc, float f, float a, int o)
{
    float noiseValue = 0.0;
    
    float localAplitude  = a;
    float localFrecuency = f;

    for( int index = 0; index < o; index++ )
    {
     	       
        noiseValue += noiseInterpolation( uv, sc * localFrecuency ) * localAplitude;
    
        localAplitude   *= 0.25;
        localFrecuency  *= 3.0;
    }    

	return noiseValue * noiseValue;
}

// =====================================================================================


void main()
{
	vec2 uv = v_TexCoords;
	vec2 suv = vec2(uv.x + 5.5, uv.y + 5.5);

	float cloudType = clamp(perlinNoise(suv, perlinScale*3.0, 0.3, 0.7,10), 0.0, 1.0);
	//float cloudType_ = clamp(perlinNoise(suv.xy, perlinScale/ 8., 0.3, 0.7,10), 0.0, 1.0);

	//cloudType = cloudType*0.5 + cloudType_*0.5;

	float coverage = perlinNoise(uv, perlinScale*0.95, perlinFrequency, perlinAmplitude, 4);
	vec4 weather = vec4( clamp(coverage, 0, 1), cloudType, 0, 1);

	//imageStore (outWeatherTex, pixel, weather);
    o_Noise = weather;
}