#pragma once

#include <iostream>

#include "Block.h"
#include "Particle.h"
#include <glm/glm.hpp>
#include "FpsCamera.h"
#include "GLClasses/VertexBuffer.h"
#include "GLClasses/VertexArray.h"
#include "GLClasses/IndexBuffer.h"
#include "GLClasses/Shader.h"
#include "BlockDatabase.h"
#include "Macros.h"

namespace VoxelRT
{
	namespace ParticleSystem
	{
		struct ParticleVertex
		{
			glm::vec3 position;
			glm::vec2 texture_coords;
		};

		class ParticleRenderer
		{
		public :
			ParticleRenderer();
			void StartParticleRender();
			void RenderParticle(const Particle& particle, FPSCamera* camera);
			void EndParticleRender(FPSCamera* camera);

		private :
			GLClasses::Shader m_ParticleShader;
			std::vector<ParticleVertex> m_ParticleVertices;

			GLClasses::VertexBuffer m_VBO;
			GLClasses::IndexBuffer m_IBO;
			GLClasses::VertexArray m_VAO;
		};

		class ParticleEmitter
		{
		public : 
			ParticleEmitter();
			void EmitParticlesAt(float lifetime, int num_particles, const glm::vec3& origin, 
				const glm::vec3& extent, const glm::vec3& vel, uint8_t block);
			void OnUpdateAndRender(FPSCamera* camera, std::array<Block, WORLD_SIZE_X* WORLD_SIZE_Y* WORLD_SIZE_Z>& data);
			void CleanUpList();

		private :
			ParticleRenderer m_Renderer;
			std::vector<Particle> m_Particles;
		};
	}
}  