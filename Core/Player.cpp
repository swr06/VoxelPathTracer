#include "Player.h"

namespace VoxelRT
{
	Player::Player() : Camera(60.0f, 800.0f / 600.0f, 0.1f, 1000.0f)
	{

	}

	bool Player::OnUpdate(GLFWwindow* window, World* world, float dt)
	{
		Camera.SetSensitivity(Sensitivity);

		bool moved = false;
		float camera_speed = Speed * dt * 4.0f;

		Camera.ResetAcceleration();
		FPSCamera cam = Camera;

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		{
			// Take the cross product of the camera's right and up.
			glm::vec3 front = -glm::cross(Camera.GetRight(), Camera.GetUp());
			Camera.ApplyAcceleration(front * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		{
			glm::vec3 back = glm::cross(Camera.GetRight(), Camera.GetUp());
			Camera.ApplyAcceleration(back * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		{
			Camera.ApplyAcceleration(-(Camera.GetRight() * camera_speed));
		}

		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		{
			Camera.ApplyAcceleration(Camera.GetRight() * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		{
			Camera.ApplyAcceleration(Camera.GetUp() * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		{
			Camera.ApplyAcceleration(-(Camera.GetUp() * camera_speed));
		}

		Camera.OnUpdate();
		glm::vec3 new_pos = Camera.GetPosition();
		glm::vec3 old_pos = cam.GetPosition();

		if (new_pos != old_pos)
		{
			glm::vec3* camera_pos = (glm::vec3*)&Camera.GetPosition();

			if (TestBlockCollision(glm::vec3(new_pos.x, old_pos.y, old_pos.z), world))
			{
				camera_pos->x = old_pos.x;
				Camera.ResetVelocity();
				Camera.ResetAcceleration();
			}

			if (TestBlockCollision(glm::vec3(old_pos.x, old_pos.y, new_pos.z), world))
			{
				camera_pos->z = old_pos.z;
				Camera.ResetVelocity();
				Camera.ResetAcceleration();
			}

			if (TestBlockCollision(glm::vec3(old_pos.x, new_pos.y, old_pos.z), world))
			{
				camera_pos->y = old_pos.y;
				Camera.ResetVelocity();
				Camera.ResetAcceleration();
			}
		}

		Camera.Refresh();

		glm::ivec3 player_block = {
			(int)floor(Camera.GetPosition().x),
			(int)floor(Camera.GetPosition().y),
			(int)floor(Camera.GetPosition().z)
		};

		InWater = false;

		return moved;
	}

	static bool Test3DAABBCollision(const glm::vec3& pos_1, const glm::vec3& dim_1, const glm::vec3& pos_2, const glm::vec3& dim_2)
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

	bool Player::TestBlockCollision(const glm::vec3& position, World* world)
	{
		if (Freefly) { return false; }

		// Convert center position to top-left position
		glm::vec3 pos = glm::vec3(
			position.x - 0.375f,
			position.y - 0.96f,
			position.z - 0.375f);

		glm::ivec3 player_block = {
			(int)floor(pos.x),
			(int)floor(pos.y),
			(int)floor(pos.z)
		};

		const glm::ivec3 block_range = { 2, 2, 2 };

		for (int i = player_block.x - block_range.x; i < player_block.x + block_range.x; i++)
			for (int j = player_block.y - block_range.y; j < player_block.y + block_range.y; j++)
				for (int k = player_block.z - block_range.z; k < player_block.z + block_range.z; k++)
				{
					if (i < WORLD_SIZE_X && i >= 0 && j < WORLD_SIZE_Y && j >= 0 && k < WORLD_SIZE_X && k >= 0)
					{
						Block* block = (Block*)&world->GetBlock(i, j, k);

						if (block && block->block != 0)
						{
							if (Test3DAABBCollision(pos, glm::vec3(0.75f, 1.5f, 0.75f), glm::vec3(i, j, k), glm::vec3(1, 1, 1)))
							{
								return true;
							}
						}
					}
				}

		return false;

	}
}