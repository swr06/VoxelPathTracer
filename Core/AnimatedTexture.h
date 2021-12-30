#pragma once

#include <glad/glad.h>

#include <iostream>
#include <vector>

namespace VoxelRT {

	class AnimatedTexture {

	public :

		GLuint m_ID;

		void Create(const std::string& path, int res, int frames);
	};


}