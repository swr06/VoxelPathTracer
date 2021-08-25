#pragma once

#include <iostream>
#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <glm/glm.hpp>
#include "GLClasses/ComputeShader.h"
#include "World.h"

namespace VoxelRT
{
	namespace Volumetrics {

		void CreateVolume(World* world);
		void PropogateVolume();
		uint8_t GetLightValue(const glm::ivec3& p);
		void SetLightValue(const glm::ivec3& p, uint8_t v);
		void AddLightToVolume(const glm::ivec3& p);
		void Reupload();
		GLuint GetVolume();
	}
}