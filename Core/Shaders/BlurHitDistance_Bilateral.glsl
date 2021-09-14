#version 330 core

layout (location = 0) out float o_HitDistance;

in vec2 v_Texcoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_HitDist;

void main() {
	float TotalWeight = 0.0f;
	float TotalDist = 0.0f;
	float BaseHitDistance = texture(u_PositionTexture,v_Texcoords).x;
	vec2 TexelSize = 1.0f/textureSize(u_HitDist,0);
	float Scale = 1.1f;

	for (int x = -2 ; x <= 2 ; x++) {
		for (int y = -2 ; y <= 2 ; y++) {
			float DistAt = texture(u_PositionTexture,v_Texcoords+vec2(x,y)*TexelSize*Scale).x;
			float HitDistAt = texture(u_HitDist,v_Texcoords+vec2(x,y)*TexelSize*Scale).x;
			float w = 1.0f / abs(DistAt-BaseHitDistance);
			w = w*w*w*w*w;
			TotalDist += HitDistAt*w;
			TotalWeight += w;
		}
	}

	o_HitDistance = TotalDist / max(TotalWeight, 0.001f);
}