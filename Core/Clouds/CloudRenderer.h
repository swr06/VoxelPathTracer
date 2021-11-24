#pragma once

#include <iostream>
#include "CloudFBO.h"
#include "NoiseRenderer.h"
#include "NoiseTexture3D.h"
#include "../FpsCamera.h"
#include "../GLClasses/Framebuffer.h"

namespace Clouds
{
	namespace CloudRenderer
	{
		void Initialize();
		GLuint Update(VoxelRT::FPSCamera& MainCamera,  
					const glm::mat4& PrevProjection,
				    const glm::mat4& PrevView, 
					const glm::vec3& CurrentPosition, 
					const glm::vec3& PrevPosition,
					GLClasses::VertexArray& VAO, const glm::vec3& SunDirection, GLuint, 
					int, int, int, GLuint atmosphere, GLuint pos_tex,
					glm::vec3 PreviousPosition, GLuint pos_tex_prev, 
					glm::vec2 modifiers, bool Clamp, glm::vec3 DetailParams, float, bool, float cirrusstrength,
					float CirrusScale, glm::ivec3 StepCounts, bool CHECKER_STEP_COUNT, float SunVisibility, 
					float CloudDetailFBMPower, bool lodlighting, bool CloudForceSupersample, float, bool);
		void RecompileShaders();
		void SetCoverage(float v);
		void SetBayer(bool v);
		void SetDetailContribution(float v);
		float GetBoxSize();
		void SetQuality(bool v);
		void SetResolution(float v);
	}
}