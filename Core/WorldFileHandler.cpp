#include "WorldFileHandler.h"

#include <sstream>
#include <filesystem>

namespace VoxelRT
{
	bool SaveWorld(World* world, const std::string& world_name)
	{
		if (!std::filesystem::exists("Saves/"))
		{
			// Create the new folder
			std::filesystem::create_directories("Saves/");
		}

		std::string world_file_name = "Saves/" + (world_name);
		FILE* world_file = NULL;

		world_file = fopen(world_file_name.c_str(), "wb+");

		if (world_file)
		{
			std::cout << "\n\n" << "SUCCESSFULLY SAVED WORLD" << "\n\n";
			fwrite(world->m_ChunkData.data(), sizeof(Block), WORLD_SIZE_X * WORLD_SIZE_Y * WORLD_SIZE_Z, world_file);
			fclose(world_file);

			return true;
		}

		else
		{
			std::cout << "\n\n" << "COULD NOT SAVE WORLD" << "\n\n";
		}

		return false;
	}

	bool LoadWorld(World* world, const std::string& world_name)
	{
		std::string world_file_name = "Saves/" + (world_name);
		FILE* world_file = NULL;

		world_file = fopen(world_file_name.c_str(), "rb");

		if (world_file)
		{
			std::cout << "\n\n" << "SUCCESSFULLY OPENED WORLD FILE" << "\n\n";
			fread(world->m_ChunkData.data(), sizeof(Block), WORLD_SIZE_X * WORLD_SIZE_Y * WORLD_SIZE_Z, world_file);
			fclose(world_file);
			std::cout << "\n\n" << "SUCCESSFULLY PARSED AND READ WORLD FILE" << "\n\n";

			return true;
		}

		else
		{
			std::cout << "\n\n" << "COULD NOT LOAD WORLD" << "\n\n";
		}

		return false;
	}

	bool FilenameValid(const std::string& name)
	{
		FILE* file = fopen(name.c_str(), "w+");

		if (!file)
		{
			return false;
		}

		fclose(file);
		std::filesystem::remove(name);

		return true;
	}
}