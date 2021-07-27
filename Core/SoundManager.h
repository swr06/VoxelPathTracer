#pragma once

#include <irrKlang.h>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "BlockDatabase.h"

namespace VoxelRT
{
	namespace SoundManager {

		void InitializeSoundManager();
		void UpdatePosition(const glm::vec3& Front, const glm::vec3& Position, const glm::vec3& Up);
		void PlaySound(const std::string& s, const glm::vec3& p, float d, float v, bool pause);
		void PlayBlockSound(uint8_t block, const glm::vec3& p, bool type);
		void Destroy();
	}
}