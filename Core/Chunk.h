#pragma once

#include <iostream>
#include <string>
#include <array>

#include "Block.h"
#include "Texture3D.h"

#define CHUNK_SIZE_X 64
#define CHUNK_SIZE_Y 64
#define CHUNK_SIZE_Z 64

namespace VoxelRT
{
	class Chunk
	{
	public :

		Chunk()
		{
			memset(&m_ChunkData, 0, CHUNK_SIZE_X * CHUNK_SIZE_Y * CHUNK_SIZE_Z);
		}

		const Block& GetBlock(uint16_t x, uint16_t y, uint16_t z)
		{
			return m_ChunkData[x + y * CHUNK_SIZE_X + z * CHUNK_SIZE_X * CHUNK_SIZE_Y];
		}

		void SetBlock(uint16_t x, uint16_t y, uint16_t z, Block block)
		{
			m_ChunkData[x + y * CHUNK_SIZE_X + z * CHUNK_SIZE_X * CHUNK_SIZE_Y] = block;
		}

		void Buffer()
		{
			m_DataTexture.CreateTexture(CHUNK_SIZE_X, CHUNK_SIZE_Y, CHUNK_SIZE_Z, m_ChunkData.data());
		}

		std::array<Block, CHUNK_SIZE_X * CHUNK_SIZE_Y * CHUNK_SIZE_Z> m_ChunkData;
		Texture3D m_DataTexture;
	};
}