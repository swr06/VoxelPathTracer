#pragma once

#include <iostream>
#include <fstream>
#include "World.h"

namespace VoxelRT
{
	bool SaveWorld(World* world, const std::string& world_name);
	bool LoadWorld(World* world, const std::string& world_name);
	bool FilenameValid(const std::string& name);
}