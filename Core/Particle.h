#pragma once

#include <iostream>
#include <chrono>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <array>
#include "Macros.h"

namespace VoxelRT
{
	namespace ParticleSystem
	{
		constexpr float gravity = 2.0;

		enum class ParticleDirection
		{
			left = 0,
			right
		};

		struct Particle
		{
			Particle(const glm::vec3& blockpos, const glm::vec3& position, const glm::vec3& velocity, const float lifetime, const float scale, ParticleDirection dir) : m_InitialPosition(position)
			{
				m_ElapsedTime = 0.0f;
				m_Rotation = 0.0f;
				m_Position = position;
				m_Velocity = velocity;
				m_Lifetime = lifetime;
				m_Scale = scale;
				m_IsAlive = true;
				m_Dir = dir;
				m_HasCollided = false;
				m_InitialPosition = blockpos;
			}

			bool TestParticleCollision(const glm::vec3& pos, std::array<Block, WORLD_SIZE_X* WORLD_SIZE_Y* WORLD_SIZE_Z>& data)
			{
				glm::ivec3 SamplePos = glm::ivec3(floor(pos.x), floor(pos.y), floor(pos.z));
				glm::ivec3 SamplePos1 = glm::ivec3(floor(pos.x), ceil(pos.y), floor(pos.z));

				if (SamplePos.x > 0 && SamplePos.x < WORLD_SIZE_X &&
					SamplePos.y > 0 && SamplePos.y < WORLD_SIZE_Y &&
					SamplePos.z > 0 && SamplePos.z < WORLD_SIZE_Z &&
					SamplePos1.x > 0 && SamplePos1.x < WORLD_SIZE_X &&
					SamplePos1.y > 0 && SamplePos1.y < WORLD_SIZE_Y &&
					SamplePos1.z > 0 && SamplePos1.z < WORLD_SIZE_Z)
				{
					uint8_t block1 = data.at(SamplePos.x + SamplePos.y * WORLD_SIZE_X + SamplePos.z * WORLD_SIZE_X * WORLD_SIZE_Y).block;
					uint8_t block2 = data.at(SamplePos1.x + SamplePos1.y * WORLD_SIZE_X + SamplePos1.z * WORLD_SIZE_X * WORLD_SIZE_Y).block;
					
					if (block1 == 0 && block2 == 0)
					{
						return false;
					}

					else
					{
						return true;
					}
				}
			}


			void OnUpdate(std::array<Block, WORLD_SIZE_X* WORLD_SIZE_Y* WORLD_SIZE_Z>& data)
			{
				glm::vec3 pos_before = m_Position;

				// Update delta every frame
				float delta = 0.1f;

				// Update the particle
				m_Velocity.y -= gravity * delta; 
				glm::vec3 change = m_Velocity;
				change *= delta;

				m_Position.y += change.y;

				float multiplier = 1.0f;

				if (TestParticleCollision(m_Position, data))
				{
					m_Position = pos_before;
					float t = 1.0f - (1.0f / (m_Lifetime - m_ElapsedTime));
					multiplier = (1.0 - t) * 0.1f + t * 0.015f;
				}

				if (m_Dir == ParticleDirection::right)
				{
					m_Position.x += change.x * multiplier;
					m_Position.z += change.z * multiplier;
				}

				else
				{
					m_Position.x -= change.x * multiplier;
					m_Position.z -= change.z * multiplier;
				}
				
				m_ElapsedTime += delta;
				m_IsAlive = IsAlive();
			}

			const glm::vec3& GetPosition() const { return m_Position; }
			const glm::vec3& GetVelocity() const { return m_Velocity; }
			float GetLifetime() const { return m_Lifetime; }
			float GetElapsedTime() const { return m_ElapsedTime; }
			const float GetScale() const { return m_Scale; }
			bool IsAlive() const { return m_Lifetime > m_ElapsedTime; }
			bool ISDead() const { return m_Lifetime < m_ElapsedTime; }

			glm::vec3 m_Position;
			glm::vec3 m_Velocity;
			glm::vec3 m_InitialPosition;

			float m_Lifetime;	
			float m_ElapsedTime;
			float m_Rotation;
			float m_Scale;
			bool m_IsAlive;
			bool m_HasCollided = false;

			uint8_t m_BlockType;
			ParticleDirection m_Dir;
		};
	}
}