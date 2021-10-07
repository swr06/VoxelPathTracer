#version 330 core

layout (location = 2) out float o_HitDistance; //location 2 -->

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_HitDist;

void main() {
	float TotalWeight = 0.0f;
	float TotalDist = 0.0f;
	float BaseHitDistance = texture(u_PositionTexture,v_TexCoords).x;
	vec2 TexelSize = 1.0f/textureSize(u_HitDist,0);
	float Scale = 1.0f;
	const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

	for (int x = -1 ; x <= 1 ; x++) {
		for (int y = -2 ; y <= 2 ; y++) {
			float DistAt = texture(u_PositionTexture,v_TexCoords+vec2(x,y)*TexelSize*Scale).x;
			float HitDistAt = texture(u_HitDist,v_TexCoords+vec2(x,y)*TexelSize*Scale).x;
			float e = abs(DistAt - BaseHitDistance);
			if (e > 1.0f) { continue; }
			float w = AtrousWeights[abs(x)]*AtrousWeights[abs(y)];
			TotalDist += HitDistAt*w;
			TotalWeight += w;
		}
	}

	o_HitDistance = TotalDist / max(TotalWeight, 0.001f);
}