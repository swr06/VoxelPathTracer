#version 330 core

layout (location = 0) out float o_Depth;

uniform sampler2D u_DepthFullRes;

void main() {
	ivec2 Offsets[4] = ivec2[4](ivec2(-1, -1), ivec2(-1, 1), ivec2(1, -1), ivec2(1, 1));
	ivec2 FullResCoordinate = ivec2(gl_FragCoord.xy) * 2;
	float MinDepth = 1000.0f;
	float EffectiveMinDepth = 2000.0f;

	for (int s = 0 ; s < 4 ; s++) {
		float SampleDepth = texelFetch(u_DepthFullRes, FullResCoordinate + Offsets[s], 0).x;
		float EffectiveDepth = SampleDepth;
		EffectiveDepth = EffectiveDepth < 0.0f ? 995.0f : EffectiveDepth;
		
		if (EffectiveDepth < EffectiveMinDepth) {
			MinDepth = SampleDepth;
			EffectiveMinDepth = EffectiveDepth;
		}
	}

	o_Depth = MinDepth;
}