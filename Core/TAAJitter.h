#pragma once

#include <iostream>
#include <glm/glm.hpp>
#include <glad/glad.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cmath>

namespace VoxelRT {

	void GenerateJitterStuff();
	glm::mat4 GetTAAJitterMatrix(int CurrentFrame, const glm::vec2& resolution);
	glm::vec2 GetTAAJitter(int CurrentFrame, const glm::vec2& resolution);
	glm::vec2 GetTAAJitterSecondary(int CurrentFrame);
}