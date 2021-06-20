#include "ParticleSystem.h"
#include "Utils/Random.h"

namespace VoxelRT
{
	namespace ParticleSystem
	{
		ParticleRenderer::ParticleRenderer() : m_VBO(GL_ARRAY_BUFFER)
		{
			m_ParticleShader.CreateShaderProgramFromFile("Core/Shaders/ParticleVert.glsl", "Core/Shaders/ParticleFrag.glsl");
			m_ParticleShader.CompileShaders();

			m_VAO.Bind();
			m_VBO.Bind();
			m_VBO.VertexAttribPointer(0, 3, GL_FLOAT, false, sizeof(ParticleVertex), (void*)offsetof(ParticleVertex, position));
			m_VBO.VertexAttribPointer(1, 2, GL_FLOAT, false, sizeof(ParticleVertex), (void*)offsetof(ParticleVertex, texture_coords));
			m_VAO.Unbind();
		}

		void ParticleRenderer::StartParticleRender()
		{
			m_ParticleVertices.clear();
		}

		void ParticleRenderer::RenderParticle(const Particle& particle, FPSCamera* camera)
		{
			glm::mat4 view_matrix = camera->GetViewMatrix();
			glm::mat4 model_matrix = glm::translate(glm::mat4(1.0f), glm::vec3(particle.m_Position));
			ParticleVertex v1, v2, v3, v4;

			// Use the transpose of the view matrix to get rid of the rotation 
			// So that the particles always face the camera
			model_matrix[0][0] = view_matrix[0][0];
			model_matrix[0][1] = view_matrix[1][0];
			model_matrix[0][2] = view_matrix[2][0];
			model_matrix[1][0] = view_matrix[0][1];
			model_matrix[1][1] = view_matrix[1][1];
			model_matrix[1][2] = view_matrix[2][1];
			model_matrix[2][0] = view_matrix[0][2];
			model_matrix[2][1] = view_matrix[1][2];
			model_matrix[2][2] = view_matrix[2][2];

			model_matrix = model_matrix * glm::scale(glm::mat4(1.0f), glm::vec3(particle.GetScale()));

			v1.position = model_matrix * glm::vec4(0.0f, 1.0f, 1.0f, 1.0f);
			v2.position = model_matrix * glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);
			v3.position = model_matrix * glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
			v4.position = model_matrix * glm::vec4(1.0f, 0.0f, 1.0f, 1.0f);

			v1.texture_coords = glm::vec2(0.0f, 1.0f);
			v2.texture_coords = glm::vec2(0.0f, 0.0f);
			v3.texture_coords = glm::vec2(1.0f, 1.0f);
			v4.texture_coords = glm::vec2(1.0f, 0.0f);

			m_ParticleVertices.push_back(v1);
			m_ParticleVertices.push_back(v2);
			m_ParticleVertices.push_back(v3);
			m_ParticleVertices.push_back(v2);
			m_ParticleVertices.push_back(v3);
			m_ParticleVertices.push_back(v4);
		}

		void ParticleRenderer::EndParticleRender(FPSCamera* camera)
		{
			if (m_ParticleVertices.size() > 0)
			{
				m_ParticleShader.Use();
				m_ParticleShader.SetMatrix4("u_ViewProjection", camera->GetViewProjection(), 0);
				
				m_VAO.Bind();
				m_VBO.BufferData(m_ParticleVertices.size() * sizeof(ParticleVertex), &m_ParticleVertices.front(), GL_STATIC_DRAW);
				glDrawArrays(GL_TRIANGLES, 0, m_ParticleVertices.size());
				m_VAO.Unbind();
			}

			m_ParticleVertices.clear();

			glUseProgram(0);
		}

		ParticleEmitter::ParticleEmitter()
		{

		}

		void ParticleEmitter::EmitParticlesAt(float lifetime, int num_particles, const glm::vec3& origin, const glm::vec3& extent,
			const glm::vec3& vel, uint8_t block)
		{
			Random random;

			for (int i = 0; i < num_particles; i++)
			{
				// Increment the x and z by a random amount
				float ix = random.UnsignedInt(extent.x) * 0.1f;
				float iy = random.UnsignedInt(extent.y) * 0.1f;
				float iz = random.UnsignedInt(extent.z) * 0.1f;
				ParticleDirection dir = ParticleDirection::left;

				if (random.Int(4) % 2 == 0) 
				{ 
					ix = -ix;
					iz = -iz; 
				}

				if (random.Int(2) == 0)
				{
					dir = ParticleDirection::right;
				}

				Particle p(glm::vec3(origin.x + ix, origin.y + iy, origin.z + iz), vel, lifetime, 0.1f, dir);
				p.m_BlockType = block;
				m_Particles.push_back(p);
			}
		}

		void ParticleEmitter::OnUpdateAndRender(FPSCamera* camera, std::array<Block, WORLD_SIZE_X* WORLD_SIZE_Y* WORLD_SIZE_Z>& data)
		{
			m_Renderer.StartParticleRender();

			for (int i = 0; i < m_Particles.size(); i++)
			{
				if (!m_Particles[i].m_IsAlive)
				{
					continue;
				}

				m_Particles[i].OnUpdate(data);
				m_Renderer.RenderParticle(m_Particles[i], camera);
			}

			m_Renderer.EndParticleRender(camera);
		}

		void ParticleEmitter::CleanUpList()
		{
			for (int i = 0; i < m_Particles.size(); i++)
			{
				if (!m_Particles[i].m_IsAlive)
				{
					m_Particles.erase(m_Particles.begin() + i);
				}
			}
		}
	}
}
