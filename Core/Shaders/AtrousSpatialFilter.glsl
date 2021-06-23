#version 330 core

#define THRESH 1.41421f

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;

uniform bool u_Dir; // 1 -> X, 0 -> Y (Meant to be separable)
uniform vec2 u_Dimensions;
uniform int u_Step;

// Edge stopping function
bool SampleValid(in vec2 SampleCoord, in vec3 InputPosition, in vec3 InputNormal)
{
	bool InScreenSpace = SampleCoord.x > 0.0f && SampleCoord.x < 1.0f && SampleCoord.y > 0.0f && SampleCoord.y < 1.0f;
	vec4 PositionAt = texture(u_PositionTexture, SampleCoord);
	vec3 NormalAt = texture(u_NormalTexture, SampleCoord).xyz;
	return (abs(PositionAt.z - InputPosition.z) <= THRESH) 
			&& (abs(PositionAt.x - InputPosition.x) <= THRESH) 
			&& (abs(PositionAt.y - InputPosition.y) <= THRESH) 
			&& (InputNormal == NormalAt) 
			&& (PositionAt.w > 0.0f) && (InScreenSpace);
}

void main()
{
	float Weights[5] = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);
	int SAMPLES = 5; // -5 to 5

	vec4 BlurredColor = vec4(0.0f);
	vec3 BasePosition = texture(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;

	if (u_Dir)
	{
		float TotalWeight = 0.0f;
		float TexelSizeX = 1.0f / u_Dimensions.x;
		
		for (int x = -SAMPLES ; x < SAMPLES ; x++)
		{
			float CurrentAtrousWeight = Weights[abs(x)];
			int SampleX = x; // todo : use u_Step here!

			vec2 SampleCoord = vec2(v_TexCoords.x + (SampleX * TexelSizeX), v_TexCoords.y);

			if (SampleValid(SampleCoord, BasePosition, BaseNormal))
			{
				vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;
				BlurredColor += SampleColor * CurrentAtrousWeight;
				TotalWeight += CurrentAtrousWeight;
			}
		}

		BlurredColor = BlurredColor / max(TotalWeight, 0.01f);
	}

	else 
	{	
		float TotalWeight = 0.0f;
		float TexelSizeY = 1.0f / u_Dimensions.y;
		
		for (int y = -SAMPLES ; y < SAMPLES ; y++)
		{
			float CurrentAtrousWeight = Weights[abs(y)];
			int SampleY = y; // todo : use u_Step here!

			vec2 SampleCoord = vec2(v_TexCoords.x, v_TexCoords.y + (SampleY * TexelSizeY));

			if (SampleValid(SampleCoord, BasePosition, BaseNormal))
			{
				vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;
				BlurredColor += SampleColor * CurrentAtrousWeight;
				TotalWeight += CurrentAtrousWeight;
			}
		}

		BlurredColor = BlurredColor / max(TotalWeight, 0.01f);
	}

	o_SpatialResult = BlurredColor;
}