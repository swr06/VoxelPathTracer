#pragma once

#include <iostream>
#include <glad/glad.h>

#include <random>
#include <numeric>
#include <vector>
#include <glm/glm.hpp>

namespace VoxelRT
{
	namespace SSAORenderer
	{
		GLuint GenerateSSAOKernelTexture();
		GLuint GenerateSSAONoiseTexture();


	}
}