#include "World.h"

#include "BlockDatabase.h"

void VoxelRT::World::ChangeCurrentlyHeldBlock()
{
	m_CurrentlyHeldBlock++;

	if (m_CurrentlyHeldBlock >= BlockDatabase::GetNumberOfBlocksInDatabase() || m_CurrentlyHeldBlock >= 127)
	{
		m_CurrentlyHeldBlock = 1;
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

void VoxelRT::World::Raycast(bool place, const glm::vec3& pos, const glm::vec3& dir)
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

				if (place)
				{
					position = position + normal;
				}

				position = glm::floor(position);

				if ((int)floor(position.x) >= WORLD_SIZE_X || (int)floor(position.y) >= WORLD_SIZE_Y || (int)floor(position.z) >= WORLD_SIZE_Z ||
					(int)floor(position.x) <= 0 || (int)floor(position.y) <= 0 || (int)floor(position.z) <= 0)
				{ 
					return; 
				}

				if (place)
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
					}
				}

				else
				{
					SetBlock((int)position.x, (int)position.y, (int)position.z, { 0 });

					if (m_Buffered)
					{
						uint8_t editblock = 0;

						glBindTexture(GL_TEXTURE_3D, m_DataTexture.GetTextureID());
						glTexSubImage3D(GL_TEXTURE_3D, 0, (int)position.x, (int)position.y, (int)position.z, 1, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &editblock);
					}
				}

				glBindTexture(GL_TEXTURE_3D, 0);

				return;
			}
		}
	}
}
