#pragma once

#include <iostream>
#include <string>
#include <array>
#include <glm/glm.hpp>

#include <algorithm>

#include "Block.h"
#include "Texture3D.h"
#include "Macros.h"

#include "GLClasses/ComputeShader.h"
#include "ParticleSystem.h"

#include "SoundManager.h"

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

		const Block& GetBlock(const glm::ivec3& p)
		{
			return m_WorldData[p.x + p.y * WORLD_SIZE_X + p.z * WORLD_SIZE_X * WORLD_SIZE_Y];
		}

		void SetBlock(const glm::ivec3& p, Block block)
		{
			m_WorldData[p.x + p.y * WORLD_SIZE_X + p.z * WORLD_SIZE_X * WORLD_SIZE_Y] = block;
		}

		void SetBlock(const glm::ivec3& p, uint8_t b)
		{
			Block block = { b };
			m_WorldData[p.x + p.y * WORLD_SIZE_X + p.z * WORLD_SIZE_X * WORLD_SIZE_Y] = block;
		}


		void InitializeLightList();
		void RebufferLightList();

		void InsertToLightList(const glm::vec3& x) {
			m_LightPositions.push_back(glm::vec4(x, 0.0f));
		}

		void RemoveFromLightList(const glm::vec3& x) {
			
			std::vector<glm::vec4>::iterator position = std::find(m_LightPositions.begin(), m_LightPositions.end(), glm::vec4(x, 0.0f));
			if (position != m_LightPositions.end()) // == myVector.end() means the element was not found
				m_LightPositions.erase(position);
		}

		void Buffer()
		{
			m_DataTexture.CreateTexture(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z, m_WorldData.data());
			m_Buffered = true;
		}

		void InitializeDistanceGenerator();
		void GenerateDistanceField();

		void ChangeCurrentlyHeldBlock(bool x);

		void Raycast(uint8_t op, const glm::vec3& pos, const glm::vec3& dir, const glm::vec3& acceleration, bool is_falling, float dt);
		void Update(FPSCamera* cam) {};
		void UpdateParticles(FPSCamera* cam, GLuint, GLuint, GLuint, GLuint, const glm::vec3& sdir, const glm::vec3& player_pos, const glm::vec2& dims, float dt);

		std::array<Block, WORLD_SIZE_X * WORLD_SIZE_Y * WORLD_SIZE_Z> m_WorldData;
		Texture3D m_DataTexture;

		std::string m_Name = "";

		uint8_t GetCurrentBlock() const noexcept { return m_CurrentlyHeldBlock; }

		Texture3D m_DistanceFieldTexture;
		ParticleSystem::ParticleEmitter m_ParticleEmitter;

		std::vector<glm::vec4> m_LightPositions;
		GLuint m_LightPositionSSBO=0;

	private :
		bool m_Buffered = false;
		uint8_t m_CurrentlyHeldBlock = 1;

		GLClasses::ComputeShader m_DistanceShaderX;
		GLClasses::ComputeShader m_DistanceShaderY;
		GLClasses::ComputeShader m_DistanceShaderZ;
	};
}