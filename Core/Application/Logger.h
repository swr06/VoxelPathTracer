#pragma once

#include <iostream>

namespace VoxelRT
{
	namespace Logger
	{
		void Log(const std::string& txt);
		void LogToFile(const std::string& txt);
	}
}