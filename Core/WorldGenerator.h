#pragma once

#include <FastNoise.h>
#include <glm/glm.hpp>
#include "World.h"

namespace VoxelRT
{
	void GenerateWorld(World* world, bool gen_type);
}