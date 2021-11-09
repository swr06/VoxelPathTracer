#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <glm/glm.hpp>
#include <queue>

#include "LightChunkMacros.h"


namespace VoxelRT {

	struct LightChunk {
		std::vector<glm::vec4> LightList;
	};
}