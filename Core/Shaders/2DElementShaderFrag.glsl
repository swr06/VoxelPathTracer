#version 330 core

in vec2 v_TexCoord;
out vec4 o_Color;
uniform sampler2D u_Texture;

void main()
{
	o_Color = texture(u_Texture, v_TexCoord);
	const float Epsilon = 0.02f;
	if (o_Color.a < Epsilon) { discard; }
}