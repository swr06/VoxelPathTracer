#pragma once

#include <iostream>
#include <string>
#include <array>
#include <glm/glm.hpp>

#include "Block.h"
#include "Texture3D.h"
#include "Macros.h"

#include "GLClasses/ComputeShader.h"
#include "ParticleSystem.h"

namespace VoxelRT
{
	class World
	{
	public :

		World()
		{
			memset(&m_WorldData, 0, WORLD_SIZE_X * WORLD_SIZE_Y * WORLD_SIZE_Z);
			m_Buffered = false;
		}

		const Block& GetBlock(uint16_t x, uint16_t y, uint16_t z)
		{
			return m_WorldData[x + y * WORLD_SIZE_X + z * WORLD_SIZE_X * WORLD_SIZE_Y];
		}

		void SetBlock(uint16_t x, uint16_t y, uint16_t z, Block block)
		{
			m_WorldData[x + y * WORLD_SIZE_X + z * WORLD_SIZE_X * WORLD_SIZE_Y] = block;
		}

		void Buffer()
		{
			m_DataTexture.CreateTexture(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z, m_WorldData.data());
			m_Buffered = true;
		}

		void InitializeDistanceGenerator();
		void GenerateDistanceField();

		void ChangeCurrentlyHeldBlock(bool x);

		void Raycast(bool place, const glm::vec3& pos, const glm::vec3& dir);
		void Update(FPSCamera* cam);

		std::array<Block, WORLD_SIZE_X * WORLD_SIZE_Y * WORLD_SIZE_Z> m_WorldData;
		Texture3D m_DataTexture;

		std::string m_Name = "";

		uint8_t GetCurrentBlock() const noexcept { return m_CurrentlyHeldBlock; }

		Texture3D m_DistanceFieldTexture;

	private :
		bool m_Buffered = false;
		uint8_t m_CurrentlyHeldBlock = 1;

		GLClasses::ComputeShader m_DistanceShaderX;
		GLClasses::ComputeShader m_DistanceShaderY;
		GLClasses::ComputeShader m_DistanceShaderZ;
		ParticleSystem::ParticleEmitter m_ParticleEmitter;
	};
}