#version 330 core

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;

in vec2 v_TexCoords;

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

void main()
{
	vec2 TexSize = textureSize(u_Texture, 0);
	vec2 TexelSize = 1.0f / TexSize;

	float TotalLuminance = 0.0f;

	for (int i = 0 ; i < TexSize.x ; i++)
	{
		for (int j = 0 ; j < TexSize.y ; j++)
		{
			TotalLuminance += (GetLuminance(texture(u_Texture, vec2(i, j) * TexelSize).rgb) + 1e-5);
		}
	}
	
	o_Color = vec3(TotalLuminance / floor(TexSize.x * TexSize.y));
}