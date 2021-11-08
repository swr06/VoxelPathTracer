#pragma once

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>

#include "BlockDatabaseParser.h"
#include "Block.h"
#include "GLClasses/TextureArray.h"

namespace VoxelRT
{
	namespace BlockDatabase
	{
		typedef uint8_t BlockIDType;

		enum BlockFaceType
		{
			Front = 0,
			Back,
			Top,
			Bottom,
			Left,
			Right
		};

		void Initialize();
		uint8_t GetBlockID(const std::string& block_name);
		int GetBlockTexture(const std::string& block_name, const BlockFaceType type);
		int GetBlockTexture(BlockIDType block_id, const BlockFaceType type);
		int GetBlockNormalTexture(const std::string& block_name, const BlockFaceType type);
		int GetBlockNormalTexture(BlockIDType block_id, const BlockFaceType type);
		int GetBlockPBRTexture(const std::string& block_name, const BlockFaceType type);
		int GetBlockPBRTexture(BlockIDType block_id, const BlockFaceType type);
		int GetBlockEmissiveTexture(BlockIDType block_id);
		std::string GetBlockName(BlockIDType block_id);
		std::string GetStepSound(BlockIDType block_id);
		std::string GetModifySound(BlockIDType block_id);

		bool IsBlockTransparent(BlockIDType block_id);

		int GetNumberOfBlocksInDatabase();

		GLuint GetTextureArray();
		GLuint GetNormalTextureArray();
		GLuint GetPBRTextureArray();
		GLuint GetEmissiveTextureArray();

		bool HasEmissiveTexture(BlockIDType block_id);
	}
}