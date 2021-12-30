#pragma once

#include <iostream>
#include <filesystem>
#include <vector>
#include <array>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>

#include <glm/glm.hpp>

#include "../Macros.h"
#include "../../Dependencies/enkiMI/enkimi.h"

#include "../BlockDatabase.h"

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

namespace VoxelRT
{
	namespace MCWorldImporter {
		void ImportWorld(const std::string& Filepath, void* Output);
	}
}