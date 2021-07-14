#include "World.h"

#include "BlockDatabase.h"

void VoxelRT::World::InitializeDistanceGenerator()
{
	int work_grp_cnt[3];

	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &work_grp_cnt[0]);
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &work_grp_cnt[1]);
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &work_grp_cnt[2]);

	printf("max global (total) work group counts x:%i y:%i z:%i\n",
		work_grp_cnt[0], work_grp_cnt[1], work_grp_cnt[2]);

	m_DistanceShaderX.CreateComputeShader("Core/Shaders/ManhattanDistanceX.comp");
	m_DistanceShaderX.Compile();
	m_DistanceShaderY.CreateComputeShader("Core/Shaders/ManhattanDistanceY.comp");
	m_DistanceShaderY.Compile();
	m_DistanceShaderZ.CreateComputeShader("Core/Shaders/ManhattanDistanceZ.comp");
	m_DistanceShaderZ.Compile();

	m_DistanceFieldTexture.CreateTexture(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z, nullptr);
}

void VoxelRT::World::GenerateDistanceField()
{
	std::cout << "\nGenerating Distance Field!\n";

	const int GROUP_SIZE = 32;

	// X PASS

	m_DistanceShaderX.Use();
	m_DistanceShaderX.SetVector3f("u_Dimensions", glm::vec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z));
	m_DistanceShaderX.SetInteger("u_BlockData", 1);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_3D, m_DataTexture.GetTextureID());

	glBindImageTexture(0, m_DistanceFieldTexture.GetTextureID(), 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8);
	glDispatchCompute(1, WORLD_SIZE_Y / GROUP_SIZE, WORLD_SIZE_Z / GROUP_SIZE);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	glUseProgram(0);


	// Y PASS

	m_DistanceShaderY.Use();

	glBindImageTexture(0, m_DistanceFieldTexture.GetTextureID(), 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8);
	glDispatchCompute(WORLD_SIZE_X / GROUP_SIZE, 1, WORLD_SIZE_Z / GROUP_SIZE);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	glUseProgram(0);


	// Z PASS

	m_DistanceShaderZ.Use();

	glBindImageTexture(0, m_DistanceFieldTexture.GetTextureID(), 0, GL_TRUE, 0, GL_READ_WRITE, GL_R8);
	glDispatchCompute(WORLD_SIZE_X / GROUP_SIZE, WORLD_SIZE_Y / GROUP_SIZE, 1);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	glUseProgram(0);

	std::cout << "\nFinished Generating Distance Field!\n";
}

void VoxelRT::World::ChangeCurrentlyHeldBlock(bool x)
{
	if (x)
	{
		m_CurrentlyHeldBlock++;

		if (m_CurrentlyHeldBlock >= BlockDatabase::GetNumberOfBlocksInDatabase() || m_CurrentlyHeldBlock >= 127)
		{
			m_CurrentlyHeldBlock = 1;
		}
	}

	else
	{
		if (m_CurrentlyHeldBlock > 0)
		{
			m_CurrentlyHeldBlock--;
		}

		if (m_CurrentlyHeldBlock <= 0)
		{
			m_CurrentlyHeldBlock = BlockDatabase::GetNumberOfBlocksInDatabase() - 1;
		}
	}
}

bool TestAABB3DCollision(const glm::vec3& pos_1, const glm::vec3& dim_1, const glm::vec3& pos_2, const glm::vec3& dim_2)
{
	if (pos_1.x < pos_2.x + dim_2.x &&
		pos_1.x + dim_1.x > pos_2.x &&
		pos_1.y < pos_2.y + dim_2.y &&
		pos_1.y + dim_1.y > pos_2.y &&
		pos_1.z < pos_2.z + dim_2.z &&
		pos_1.z + dim_1.z > pos_2.z)
	{
		return true;
	}

	return false;
}

bool TestRayPlayerCollision(const glm::vec3& ray_block, const glm::vec3& player_pos)
{
	glm::vec3 pos = player_pos; 

	if (TestAABB3DCollision(pos, glm::vec3(0.75f, 1.5f, 0.75f), ray_block, glm::vec3(1.2f, 1.2f, 1.2f)))
	{
		return true;
	}

	return false;

}

void VoxelRT::World::Raycast(uint8_t op, const glm::vec3& pos, const glm::vec3& dir)
{
	glm::vec3 position = pos;
	const glm::vec3& direction = dir;
	int max = 50; // block reach

	glm::vec3 sign;

	for (int i = 0; i < 3; ++i)
		sign[i] = direction[i] > 0;

	for (int i = 0; i < max; ++i)
	{
		glm::vec3 tvec = (floor(position + sign) - position) / direction;
		float t = std::min(tvec.x, std::min(tvec.y, tvec.z));

		position += direction * (t + 0.001f);

		if (!((int)floor(position.x) >= WORLD_SIZE_X || (int)floor(position.y) >= WORLD_SIZE_Y || (int)floor(position.z) >= WORLD_SIZE_Z ||
			(int)floor(position.x) <= 0 || (int)floor(position.y) <= 0 || (int)floor(position.z) <= 0 ))
		{
			Block ray_block = GetBlock((int)position.x, (int)position.y, (int)position.z);

			if (ray_block.block != 0)
			{
				glm::vec3 normal;

				for (int j = 0; j < 3; ++j)
				{
					normal[j] = (t == tvec[j]);

					if (sign[j])
					{
						normal[j] = -normal[j];
					}
				}

				if (op == 1)
				{
					position = position + normal;
				}

				position = glm::floor(position);

				if ((int)floor(position.x) >= WORLD_SIZE_X || (int)floor(position.y) >= WORLD_SIZE_Y || (int)floor(position.z) >= WORLD_SIZE_Z ||
					(int)floor(position.x) <= 0 || (int)floor(position.y) <= 0 || (int)floor(position.z) <= 0)
				{ 
					return; 
				}

				if (op == 1)
				{
					if (TestRayPlayerCollision(glm::vec3((position.x), (position.y), (position.z)), pos))
					{
						return;
					}

					uint8_t editblock = m_CurrentlyHeldBlock;

					SetBlock((int)position.x, (int)position.y, (int)position.z, { editblock });
					
					if (m_Buffered)
					{
					
						glBindTexture(GL_TEXTURE_3D, m_DataTexture.GetTextureID());
						glTexSubImage3D(GL_TEXTURE_3D, 0, (int)position.x, (int)position.y, (int)position.z, 1, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &editblock);
						GenerateDistanceField();
					}
				}

				else if (op == 0)
				{
					float y1 = rand() % 2 == 0 ? 0.05f : 0.075f;
					float y2 = y1 == 0.05f ? 0.1f : 0.113f;
					glm::vec3 x = rand() % 2 == 0 ? glm::vec3(y1, 1, y2) : glm::vec3(y2, 1, y1);
					glm::vec3 particle_pos;
					particle_pos.x = floor(position.x);
					particle_pos.y = position.y + 0.4f;
					particle_pos.z = floor(position.z);
					particle_pos.x += 0.5f;
					particle_pos.z += 0.5f;
					m_ParticleEmitter.EmitParticlesAt(glm::vec3(floor(position.x), floor(position.y), floor(position.z)), 2.75f, 40,
					particle_pos, glm::vec3(5, 5, 5), x, GetBlock((int)position.x, (int)position.y, (int)position.z).block);

					SetBlock((int)position.x, (int)position.y, (int)position.z, { 0 });

					if (m_Buffered)
					{
						uint8_t editblock = 0;

						glBindTexture(GL_TEXTURE_3D, m_DataTexture.GetTextureID());
						glTexSubImage3D(GL_TEXTURE_3D, 0, (int)position.x, (int)position.y, (int)position.z, 1, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &editblock);
						GenerateDistanceField();
					}
				}

				else if (op == 2)
				{
					uint8_t block = GetBlock((int)position.x, (int)position.y, (int)position.z).block;

					if (block > 0)
					{
						m_CurrentlyHeldBlock = block;
					}

					return;
				}

				glBindTexture(GL_TEXTURE_3D, 0);

				return;
			}
		}
	}
}

void VoxelRT::World::UpdateParticles(FPSCamera* cam, GLuint pos_texture, GLuint shadow_tex, GLuint diff, GLuint diff2, const glm::vec3& sdir, const glm::vec3& player_pos, const glm::vec2& dims, float dt)
{
	m_ParticleEmitter.OnUpdateAndRender(cam, m_WorldData, pos_texture, shadow_tex, diff, diff2, sdir, player_pos, dims, dt);
}
