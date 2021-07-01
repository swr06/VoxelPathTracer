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
					GLClasses::VertexArray& VAO, const glm::vec3& SunDirection, GLuint, int, int, int, GLuint atmosphere, glm::vec3 PreviousPosition);
		void RecompileShaders();
		void SetChecker(bool v);
		void SetCoverage(float v);
		void SetBayer(bool v);
		void SetDetailContribution(float v);
	}
}