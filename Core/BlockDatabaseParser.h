#pragma once

#include <iostream>
#include <map>
#include <unordered_map>
#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include <sstream>

#include "Application/Logger.h"

namespace VoxelRT
{
	namespace BlockDatabaseParser
	{
		struct BlockTexture
		{
			std::string front = "";
			std::string back = "";
			std::string left = "";
			std::string right = "";
			std::string top = "";
			std::string bottom = "";
		};

		struct ParsedBlockData
		{
			BlockTexture NormalMap;
			BlockTexture PBRMap;
			BlockTexture AlbedoMap;
			std::string EmissiveMap = "";
			std::string BlockName = "";
			uint8_t ID = 0;
			bool transparent = false;
			bool sss = false;
			std::string snd_step = "";
			std::string snd_modify = "";
		};

		void Parse(const std::string& path);
	}
}