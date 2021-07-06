#version 330 core
#define THRESH 1.41421f // root 2

layout (location = 0) out vec4 o_SpatialResult;

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_InputTexture;
uniform vec2 u_Dimensions;

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

vec4 BlurVertical(int x, in vec3 BasePosition, in vec3 BaseNormal) 
{
	const float[5] Weights = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);
	vec2 TexelSize = 1.0f / u_Dimensions;
	
	vec4 Total = vec4(0.0f);
	float TotalWeight = 0.0f;

	for (int y = -2; y <= 2; y++)
	{
		float CurrentWeight = Weights[y + 2];
		vec2 SampleCoord = v_TexCoords + (vec2(x, y) * TexelSize);

		if (SampleValid(SampleCoord, BasePosition, BaseNormal))
		{
			vec4 SampleColor = texture(u_InputTexture, SampleCoord).xyzw;
			Total += SampleColor * CurrentWeight;
			TotalWeight += CurrentWeight;
		}
	}

	Total = Total / max(TotalWeight, 0.01f);
	return Total;
}

void main()
{
	const float[5] Weights = float[5] (0.0625, 0.25, 0.375, 0.25, 0.0625);
	vec2 TexelSize = 1.0f / u_Dimensions;

	vec4 TotalColor = vec4(0.0f);
	float TotalWeight = 0.0f; 

	vec3 BasePosition = texture(u_PositionTexture, v_TexCoords).xyz;
	vec3 BaseNormal = texture(u_NormalTexture, v_TexCoords).xyz;

	for (int x = -2 ; x <= 2 ; x++)
	{
		float CurrentWeight = Weights[x + 2];
		TotalColor += BlurVertical(x, BasePosition, BaseNormal) * CurrentWeight;
		TotalWeight += CurrentWeight;
	}

	o_SpatialResult = TotalColor / max(TotalWeight, 0.01f);
}