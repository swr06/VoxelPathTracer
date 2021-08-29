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

		void CreateVolume(World* world, GLuint SSBO_Blockdata, GLuint AlbedoArray);
		void PropogateVolume();
		uint8_t GetLightValue(const glm::ivec3& p);
		uint8_t GetBlockTypeLightValue(const glm::ivec3& p);
		void SetLightValue(const glm::ivec3& p, uint8_t v, uint8_t block);
		void UploadLight(const glm::ivec3& p, uint8_t v, uint8_t block, bool should_bind);
		void AddLightToVolume(const glm::ivec3& p, uint8_t block);
		void Reupload();
		GLuint GetVolume();
		GLuint GetAverageColorSSBO();
	}
}