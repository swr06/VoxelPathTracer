#pragma once

#include <glad/glad.h>

#include <iostream>
#include <vector>
#include <array>

namespace VoxelRT
{
	class BlueNoiseDataSSBO
	{
	public :
		BlueNoiseDataSSBO();
		~BlueNoiseDataSSBO();

		GLuint m_SSBO;
	};

}