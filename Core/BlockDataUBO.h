#pragma once

#include <iostream>
#include <glad/glad.h>
#include "BlockDatabase.h"
#include "ShaderManager.h"

#include "Macros.h"

#include <array>
#include <vector>

namespace VoxelRT
{
	class BlockDataUBO
	{
	public :

		void CreateBuffers();
		inline GLuint GetSSBO() const noexcept { return m_SSBO; }
		void Bind(int idx);

	private :

		GLuint m_SSBO;
	};
}