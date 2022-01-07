// Properly downsamples the gbuffer


#version 330 core

#define LINEARIZE_DEPTH 


layout (location = 0) out float o_Depth;
layout (location = 1) out float o_Normals;
layout (location = 2) out float o_BlockIDs;

// Non linear high precision depth ->
uniform sampler2D u_DepthFullRes; 

uniform sampler2D u_NormalsFullRes;
uniform sampler2D u_BlockIDs;

void main() {
	ivec2 Offsets[4] = ivec2[4](ivec2(-1, -1), ivec2(-1, 1), ivec2(1, -1), ivec2(1, 1));
	
	ivec2 FullResCoordinate = ivec2(gl_FragCoord.xy) * 2;
	float EffectiveMinDepth = 2400.0f;
	ivec2 MinOffset = ivec2(0);
	float MinDepth = 0.0f;

	for (int s = 0 ; s < 4 ; s++) {

		float SampleDepth = texelFetch(u_DepthFullRes, FullResCoordinate + Offsets[s], 0).x;

		#ifdef LINEARIZE_DEPTH
			SampleDepth = 1.0f / SampleDepth;
		#endif
		
		float EffectiveDepth = SampleDepth;
		EffectiveDepth = EffectiveDepth < 0.0f ? 600.0f : EffectiveDepth;
		
		if (EffectiveDepth < EffectiveMinDepth) {
			EffectiveMinDepth = EffectiveDepth;
			MinOffset = Offsets[s];
			MinDepth = SampleDepth;
		}
	}

	float DownsampledNormal = texelFetch(u_NormalsFullRes, FullResCoordinate + MinOffset, 0).x;
	float DownsampledID = texelFetch(u_BlockIDs, FullResCoordinate + MinOffset, 0).x;
	float DownsampledDepth = MinDepth;

	o_Depth = DownsampledDepth;
	o_Normals = DownsampledNormal;
	o_BlockIDs = DownsampledID;
}