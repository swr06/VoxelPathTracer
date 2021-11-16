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
	//glm::ivec3 Get3DIdx(int idx)
	//{
	//	int z = idx / (WORLD_SIZE_X * WORLD_SIZE_Y);
	//	idx -= (z * WORLD_SIZE_X * WORLD_SIZE_Y);
	//	int y = idx / WORLD_SIZE_X;
	//	int x = idx % WORLD_SIZE_X;
	//	return glm::ivec3(x, y, z);
	//}

	inline static int Get1DIndexForLightChunk(int x, int y, int z)
	{
		return ((z * 24 * 8) + (y * 24) + x);
	}

	class World
	{
	public :

		World()
		{
			memset(&m_WorldData, 0, WORLD_SIZE_X * WORLD_SIZE_Y * WORLD_SIZE_Z);
			memset(&LightChunkOffsets[0], -1, 24 * 8 * 24 * sizeof(int) * 2);
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
		//void RebufferLightList();

		void InsertToLightList(const glm::vec3& x) {

			//m_LightPositions.push_back(glm::vec4(x, 0.0f));

			// Place light inside chunk data and update offsets
			int Cx = floor((float)x.x / (float)16); // Cx 
			int Cy = floor((float)x.y / (float)16);
			int Cz = floor((float)x.z / (float)16);
			int CurrentIDX = Get1DIndexForLightChunk(Cx, Cy, Cz); // Get chunk index
			auto& StoredOffset = LightChunkOffsets[CurrentIDX]; // get the offset stored 

			if (StoredOffset.x > 0) // if there were previously placed lights here, insert it at the stored offset
			{
				LightChunkData.insert(LightChunkData.begin() + StoredOffset.x, glm::vec4(x, 0.0f));
				StoredOffset.y = glm::max(StoredOffset.y + 1, 1); // increment size

				// Update the offset for the other light chunks
				for (int x = 0; x < LightChunkOffsets.size(); x++) 
				{
					if (x == CurrentIDX) {
						continue;
					}

					if (LightChunkOffsets[x].x >= StoredOffset.x) {
						LightChunkOffsets[x].x++; // Shift by one
					}
				}
			}

			else {
				// if this is the first light placed in the chunk, push it back to the vector and store the offset
				StoredOffset.x = LightChunkData.size(); // down perhaps?
				LightChunkData.push_back(glm::vec4(x, 0.0f));
				StoredOffset.y = 1;
			}

			RebufferLightChunks();
		}

		void RemoveFromLightList(const glm::vec3& x) {
			
			//std::vector<glm::vec4>::iterator position = std::find(m_LightPositions.begin(), m_LightPositions.end(), glm::vec4(x, 0.0f));
			//if (position != m_LightPositions.end()) // == myVector.end() means the element was not found
			//{
			//	m_LightPositions.erase(position);
			//}

			// Place light inside chunk data and update offsets
			int Cx = floor((float)x.x / (float)16); // Cx 
			int Cy = floor((float)x.y / (float)16);
			int Cz = floor((float)x.z / (float)16);
			int CurrentIDX = Get1DIndexForLightChunk(Cx, Cy, Cz); // Get chunk index
			auto& StoredOffset = LightChunkOffsets[CurrentIDX]; // get the offset stored 
			int IterationCount = 0;

			if (StoredOffset.x > 0) // if there were previously placed lights here, remove it at the stored offset
			{
				bool StopIterating = false;

				while (StopIterating == false) {

					if (IterationCount > 2) {
						break; 
					}

					IterationCount++;

					std::vector<glm::vec4>::iterator ChunkIter = std::find(LightChunkData.begin() + StoredOffset.x, LightChunkData.begin() + StoredOffset.x + StoredOffset.y, glm::vec4(x, 0.0f));
					if (ChunkIter != LightChunkData.end()) // The element was found
					{
						LightChunkData.erase(ChunkIter); // erase it
						StoredOffset.y = glm::max(StoredOffset.y - 1, 0);

						for (int x = 0; x < LightChunkOffsets.size(); x++)
						{
							if (x == CurrentIDX) {
								continue;
							}

							if (LightChunkOffsets[x].x >= StoredOffset.x) {
								LightChunkOffsets[x].x--;
								LightChunkOffsets[x].x = glm::max(LightChunkOffsets[x].x, 0);
							}
						}
					}

					else {
						StopIterating = true;
					}
				}
			}

			RebufferLightChunks();
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

		//std::vector<glm::vec4> m_LightPositions;

		void RepropogateLPV_();
		void RebufferLightChunks();

		// Each "chunk" stores a list of lights that can be sampled using MIS
		std::vector<glm::vec4> LightChunkData;
		std::array<glm::ivec2, 24 * 8 * 24> LightChunkOffsets; // 16x16x16 chunks

		GLuint LightChunkDataSSBO = 0;
		GLuint LightChunkOffsetSSBO = 0;

	private :
		bool m_Buffered = false;
		uint8_t m_CurrentlyHeldBlock = 1;

		GLClasses::ComputeShader m_DistanceShaderX;
		GLClasses::ComputeShader m_DistanceShaderY;
		GLClasses::ComputeShader m_DistanceShaderZ;
		
	};
}