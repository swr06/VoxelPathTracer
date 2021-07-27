#include "BlockDatabase.h"

namespace VoxelRT
{
	extern std::unordered_map<std::string, BlockDatabaseParser::ParsedBlockData> ParsedBlockDataList;
	std::unordered_map<uint8_t, BlockDatabaseParser::ParsedBlockData> ParsedBlockDataListID;
	GLClasses::TextureArray BlockTextureArray;
	GLClasses::TextureArray BlockNormalTextureArray;
	GLClasses::TextureArray BlockPBRTextureArray;
	GLClasses::TextureArray BlockEmissiveTextureArray;
	
	void BlockDatabase::Initialize()
	{
		Logger::Log("Starting to parse database file");
		std::string database_path = "blockdb.txt";
		BlockDatabaseParser::Parse(database_path);

		for (auto& e : ParsedBlockDataList)
		{
			uint8_t id = e.second.ID;
			const BlockDatabaseParser::ParsedBlockData& data = e.second;

			ParsedBlockDataListID[id] = data;
		}

		Logger::Log("Successfully parsed database file");

		std::vector<std::string> paths;
		std::vector<std::string> normal_paths;
		std::vector<std::string> pbr_paths;
		std::vector<std::string> emissive_paths;

		std::pair<int, int> texture_resolutions = { 512,512 };

		for (auto& e : ParsedBlockDataList)
		{
			paths.push_back(e.second.AlbedoMap.front);
			paths.push_back(e.second.AlbedoMap.back);
			paths.push_back(e.second.AlbedoMap.top);
			paths.push_back(e.second.AlbedoMap.bottom);
			paths.push_back(e.second.AlbedoMap.left);
			paths.push_back(e.second.AlbedoMap.right);

			normal_paths.push_back(e.second.NormalMap.front);
			normal_paths.push_back(e.second.NormalMap.back);
			normal_paths.push_back(e.second.NormalMap.top);
			normal_paths.push_back(e.second.NormalMap.bottom);
			normal_paths.push_back(e.second.NormalMap.left);
			normal_paths.push_back(e.second.NormalMap.right);

			pbr_paths.push_back(e.second.PBRMap.front);
			pbr_paths.push_back(e.second.PBRMap.back);
			pbr_paths.push_back(e.second.PBRMap.top);
			pbr_paths.push_back(e.second.PBRMap.bottom);
			pbr_paths.push_back(e.second.PBRMap.left);
			pbr_paths.push_back(e.second.PBRMap.right);
			
			if (e.second.EmissiveMap.size() > 0)
			{
				emissive_paths.push_back(e.second.EmissiveMap);
			}
		}

		std::string res_str = "     |     RES : (" + std::to_string(texture_resolutions.first) + "," + std::to_string(texture_resolutions.second) + ")";

		std::cout << ("\n\n\n-- STARTING TO LOAD TEXTURES --\n\n\n");

		Logger::Log("Creating Albedo texture array!" + res_str);
		BlockTextureArray.CreateArray(paths, texture_resolutions, true, true);
		Logger::Log("Successfully created Albedo texture array!");

		Logger::Log("Creating Normal texture array!" + res_str);
		BlockNormalTextureArray.CreateArray(normal_paths, texture_resolutions, false, true, GL_LINEAR, true);
		Logger::Log("Successfully created Albedo texture array!");

		Logger::Log("Creating PBR texture array!" + res_str);
		BlockPBRTextureArray.CreateArray(pbr_paths, texture_resolutions, false, true, GL_LINEAR, true);
		Logger::Log("Successfully created PBR texture array!");

		Logger::Log("Creating Emissive texture array!" + res_str);
		BlockEmissiveTextureArray.CreateArray(emissive_paths, texture_resolutions, false, true, GL_LINEAR, true);
		Logger::Log("Successfully created emissive texture array!");

		std::cout << ("\n\n\n-- SUCCESSFULLY LOADED TEXTURES --\n\n\n");
	}

	uint8_t BlockDatabase::GetBlockID(const std::string& block_name)
	{
		return ParsedBlockDataList[block_name].ID;
	}

	int BlockDatabase::GetBlockTexture(const std::string& block_name, const BlockFaceType type)
	{
		if (ParsedBlockDataList.find(block_name) == ParsedBlockDataList.end())
		{
			return -1;
		}

		else
		{
			std::string pth;

			switch (type)
			{
				case BlockFaceType::Front : 
				{	
					pth = ParsedBlockDataList[block_name].AlbedoMap.front;
					break;
				}

				case BlockFaceType::Back:
				{
					pth = ParsedBlockDataList[block_name].AlbedoMap.back;
					break;
				}

				case BlockFaceType::Left:
				{
					pth = ParsedBlockDataList[block_name].AlbedoMap.left;
					break;
				}

				case BlockFaceType::Right:
				{
					pth = ParsedBlockDataList[block_name].AlbedoMap.right;
					break;
				}

				case BlockFaceType::Top:
				{
					pth = ParsedBlockDataList[block_name].AlbedoMap.top;
					break;
				}

				case BlockFaceType::Bottom:
				{
					pth = ParsedBlockDataList[block_name].AlbedoMap.bottom;
					break;
				}
			}

			return BlockTextureArray.GetTexture(pth);
		}

		return -1;
	}

	int BlockDatabase::GetBlockTexture(BlockIDType block_id, const BlockFaceType type)
	{
		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return -1;
		}

		else
		{
			std::string pth;

			switch (type)
			{
			case BlockFaceType::Front:
			{
				pth = ParsedBlockDataListID[block_id].AlbedoMap.front;
				break;
			}

			case BlockFaceType::Back:
			{
				pth = ParsedBlockDataListID[block_id].AlbedoMap.back;
				break;
			}

			case BlockFaceType::Left:
			{
				pth = ParsedBlockDataListID[block_id].AlbedoMap.left;
				break;
			}

			case BlockFaceType::Right:
			{
				pth = ParsedBlockDataListID[block_id].AlbedoMap.right;
				break;
			}

			case BlockFaceType::Top:
			{
				pth = ParsedBlockDataListID[block_id].AlbedoMap.top;
				break;
			}

			case BlockFaceType::Bottom:
			{
				pth = ParsedBlockDataListID[block_id].AlbedoMap.bottom;
				break;
			}
			}

			return BlockTextureArray.GetTexture(pth);
		}

		return -1;
	}

	// Normal
	int BlockDatabase::GetBlockNormalTexture(const std::string& block_name, const BlockFaceType type)
	{
		if (ParsedBlockDataList.find(block_name) == ParsedBlockDataList.end())
		{
			return -1;
		}

		else
		{
			std::string pth;

			switch (type)
			{
			case BlockFaceType::Front:
			{
				pth = ParsedBlockDataList[block_name].NormalMap.front;
				break;
			}

			case BlockFaceType::Back:
			{
				pth = ParsedBlockDataList[block_name].NormalMap.back;
				break;
			}

			case BlockFaceType::Left:
			{
				pth = ParsedBlockDataList[block_name].NormalMap.left;
				break;
			}

			case BlockFaceType::Right:
			{
				pth = ParsedBlockDataList[block_name].NormalMap.right;
				break;
			}

			case BlockFaceType::Top:
			{
				pth = ParsedBlockDataList[block_name].NormalMap.top;
				break;
			}

			case BlockFaceType::Bottom:
			{
				pth = ParsedBlockDataList[block_name].NormalMap.bottom;
				break;
			}
			}

			return BlockNormalTextureArray.GetTexture(pth);
		}

		return -1;
	}

	// Normal
	int BlockDatabase::GetBlockNormalTexture(BlockIDType block_id, const BlockFaceType type)
	{
		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return -1;
		}

		else
		{
			std::string pth;

			switch (type)
			{
			case BlockFaceType::Front:
			{
				pth = ParsedBlockDataListID[block_id].NormalMap.front;
				break;
			}

			case BlockFaceType::Back:
			{
				pth = ParsedBlockDataListID[block_id].NormalMap.back;
				break;
			}

			case BlockFaceType::Left:
			{
				pth = ParsedBlockDataListID[block_id].NormalMap.left;
				break;
			}

			case BlockFaceType::Right:
			{
				pth = ParsedBlockDataListID[block_id].NormalMap.right;
				break;
			}

			case BlockFaceType::Top:
			{
				pth = ParsedBlockDataListID[block_id].NormalMap.top;
				break;
			}

			case BlockFaceType::Bottom:
			{
				pth = ParsedBlockDataListID[block_id].NormalMap.bottom;
				break;
			}
			}

			return BlockNormalTextureArray.GetTexture(pth);
		}

		return -1;
	}

	// PBR

	int BlockDatabase::GetBlockPBRTexture(const std::string& block_name, const BlockFaceType type)
	{
		if (ParsedBlockDataList.find(block_name) == ParsedBlockDataList.end())
		{
			return -1;
		}

		else
		{
			std::string pth;

			switch (type)
			{
			case BlockFaceType::Front:
			{
				pth = ParsedBlockDataList[block_name].PBRMap.front;
				break;
			}

			case BlockFaceType::Back:
			{
				pth = ParsedBlockDataList[block_name].PBRMap.back;
				break;
			}

			case BlockFaceType::Left:
			{
				pth = ParsedBlockDataList[block_name].PBRMap.left;
				break;
			}

			case BlockFaceType::Right:
			{
				pth = ParsedBlockDataList[block_name].PBRMap.right;
				break;
			}

			case BlockFaceType::Top:
			{
				pth = ParsedBlockDataList[block_name].PBRMap.top;
				break;
			}

			case BlockFaceType::Bottom:
			{
				pth = ParsedBlockDataList[block_name].PBRMap.bottom;
				break;
			}
			}

			return BlockPBRTextureArray.GetTexture(pth);
		}

		return -1;
	}

	// PBR
	int BlockDatabase::GetBlockPBRTexture(BlockIDType block_id, const BlockFaceType type)
	{
		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return 0;
		}

		else
		{
			std::string pth;

			switch (type)
			{
			case BlockFaceType::Front:
			{
				pth = ParsedBlockDataListID[block_id].PBRMap.front;
				break;
			}

			case BlockFaceType::Back:
			{
				pth = ParsedBlockDataListID[block_id].PBRMap.back;
				break;
			}

			case BlockFaceType::Left:
			{
				pth = ParsedBlockDataListID[block_id].PBRMap.left;
				break;
			}

			case BlockFaceType::Right:
			{
				pth = ParsedBlockDataListID[block_id].PBRMap.right;
				break;
			}

			case BlockFaceType::Top:
			{
				pth = ParsedBlockDataListID[block_id].PBRMap.top;
				break;
			}

			case BlockFaceType::Bottom:
			{
				pth = ParsedBlockDataListID[block_id].PBRMap.bottom;
				break;
			}
			}

			return BlockPBRTextureArray.GetTexture(pth);
		}

		return -1;
	}

	int BlockDatabase::GetBlockEmissiveTexture(BlockIDType block_id)
	{
		std::string pth;

		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return -1;
		}

		pth = ParsedBlockDataListID[block_id].EmissiveMap;

		if (pth.size() > 0)
		{
			return BlockEmissiveTextureArray.GetTexture(pth);
		}

		return -1;
	}

	std::string BlockDatabase::GetBlockName(BlockIDType block_id)
	{
		if (block_id == 0) { return "???"; }

		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return "???";
		}

		return ParsedBlockDataListID.at(block_id).BlockName;
	}

	std::string BlockDatabase::GetStepSound(BlockIDType block_id)
	{
		if (block_id == 0) { return "???"; }

		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return "???";
		}

		return ParsedBlockDataListID.at(block_id).snd_step;
	}

	std::string BlockDatabase::GetModifySound(BlockIDType block_id)
	{
		if (block_id == 0) { return "???"; }

		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return "???";
		}

		return ParsedBlockDataListID.at(block_id).snd_modify;
	}

	bool BlockDatabase::IsBlockTransparent(BlockIDType block_id)
	{
		if (block_id == 0) { return true; }

		if (ParsedBlockDataListID.find(block_id) == ParsedBlockDataListID.end())
		{
			return false;
		}

		return ParsedBlockDataListID.at(block_id).transparent;
	}

	int BlockDatabase::GetNumberOfBlocksInDatabase()
	{
		return ParsedBlockDataList.size();
	}

	GLuint BlockDatabase::GetTextureArray()
	{
		return BlockTextureArray.GetTextureArray();
	}

	GLuint BlockDatabase::GetNormalTextureArray()
	{
		return BlockNormalTextureArray.GetTextureArray();
	}

	GLuint BlockDatabase::GetPBRTextureArray()
	{
		return BlockPBRTextureArray.GetTextureArray();
	}

	GLuint BlockDatabase::GetEmissiveTextureArray()
	{
		return BlockEmissiveTextureArray.GetTextureArray();
	}
}
