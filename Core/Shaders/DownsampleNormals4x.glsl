#version 330 core

layout (location = 0) out float o_NormalID;

uniform sampler2D u_DepthFullRes;
uniform sampler2D u_NormalsFullRes;

void main() {
	// since the resolution is 4x lesser, we need to sample 4 texels 
	ivec2 Offsets[4] = ivec2[4](ivec2(-1, -1), ivec2(-1, 1), ivec2(1, -1), ivec2(1, 1));
	
	ivec2 FullResCoordinate = ivec2(gl_FragCoord.xy) * 2;
	float EffectiveMinDepth = 2000.0f;
	ivec2 MinOffset = ivec2(0);

	for (int s = 0 ; s < 4 ; s++) {

		float SampleDepth = texelFetch(u_DepthFullRes, FullResCoordinate + Offsets[s], 0).x;
		float EffectiveDepth = SampleDepth;
		EffectiveDepth = EffectiveDepth < 0.0f ? 400.0f : EffectiveDepth;
		
		if (EffectiveDepth < EffectiveMinDepth) {
			EffectiveMinDepth = EffectiveDepth;
			MinOffset = Offsets[s];
		}
	}

	o_NormalID = texelFetch(u_NormalsFullRes, FullResCoordinate + MinOffset, 0).x;
}