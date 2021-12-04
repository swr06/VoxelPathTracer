#version 330 core

layout (location = 0) out float o_HitDistance;

in vec2 v_TexCoords;

uniform sampler2D u_PositionTexture;
uniform sampler2D u_HitDist;

void main() {
	
	
	vec2 TexelSize = 1.0f/textureSize(u_HitDist,0);
	float Scale = 4.0f;
	const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );


	float BaseHitDistance = texture(u_PositionTexture,v_TexCoords).x;
	float TotalWeight = 1.0f;
	float TotalDist = BaseHitDistance;
	int SamplesValid = 0;


	for (int x = -2 ; x <= 2 ; x++) {
		for (int y = -2 ; y <= 2 ; y++) {
			if (x == 0 && y == 0) { continue; }
			float DistAt = texture(u_PositionTexture,v_TexCoords+vec2(x,y)*TexelSize*Scale).x;
			float HitDistAt = texture(u_HitDist,v_TexCoords+vec2(x,y)*TexelSize*Scale).x;
			float e = abs(DistAt - BaseHitDistance);
			if (e > 1.0f || HitDistAt < 0.0f) { continue; }
			float w = AtrousWeights[abs(x)]*AtrousWeights[abs(y)];
			TotalDist += HitDistAt*w;
			TotalWeight += w;
			SamplesValid++;
		}
	}

	if (SamplesValid == 0) {
		TotalWeight = 1.0f;
	    TotalDist = BaseHitDistance;
	}

	o_HitDistance = TotalDist / max(TotalWeight, 0.001f);
}