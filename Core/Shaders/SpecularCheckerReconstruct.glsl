#version 330 core

layout (location = 0) out vec4 o_SH;
layout (location = 1) out vec2 o_CoCg;

in vec2 v_TexCoords;
uniform int u_CurrentFrame;
uniform sampler2D u_CurrentSH;
uniform sampler2D u_CoCgTexture;

vec4 TextureClamped(sampler2D t, vec2 txc) {
	return texture(t, clamp(txc, 0.0f, 1.0f));
}


void main()
{
	int CheckerStep;
	CheckerStep = u_CurrentFrame % 2;
	vec2 TexelSize = 1.0f / textureSize(u_CurrentSH, 0).xy;

	if (int(gl_FragCoord.x + gl_FragCoord.y) % 2 == CheckerStep)
	{	
		const ivec2 Offsets[4] = ivec2[](ivec2(1, 0), ivec2(0, 1), ivec2(-1, 0), ivec2(0, -1));

		vec4 Total = vec4(0.0f);
		Total += TextureClamped(u_CurrentSH, v_TexCoords + (Offsets[0] * TexelSize));
		Total += TextureClamped(u_CurrentSH, v_TexCoords + (Offsets[1] * TexelSize));
		Total += TextureClamped(u_CurrentSH, v_TexCoords + (Offsets[2] * TexelSize));
		Total += TextureClamped(u_CurrentSH, v_TexCoords + (Offsets[3] * TexelSize));
		Total /= 4.0f;


		vec2 TotalCoCg = vec2(0.0f);
		TotalCoCg += TextureClamped(u_CoCgTexture, v_TexCoords + (Offsets[0] * TexelSize)).rg;
		TotalCoCg += TextureClamped(u_CoCgTexture, v_TexCoords + (Offsets[1] * TexelSize)).rg;
		TotalCoCg += TextureClamped(u_CoCgTexture, v_TexCoords + (Offsets[2] * TexelSize)).rg;
		TotalCoCg += TextureClamped(u_CoCgTexture, v_TexCoords + (Offsets[3] * TexelSize)).rg;
		TotalCoCg /= 4.0f;

		o_SH = Total;
		o_CoCg = TotalCoCg;
	}	

	else 
	{
		vec4 Fetch = texture(u_CurrentSH, v_TexCoords);
		vec2 CoCg = texture(u_CoCgTexture, v_TexCoords).rg;

		o_SH = Fetch;
		o_CoCg = CoCg;
	}
}