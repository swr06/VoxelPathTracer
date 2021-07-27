#include "Player.h"

namespace VoxelRT
{
	Player::Player() : Camera(60.0f, 800.0f / 600.0f, 0.1f, 1000.0f), m_AABB(glm::vec3(0.3f, 1.0f, 0.3f))
	{
		m_Acceleration = glm::vec3(0.0f);
		m_Velocity = glm::vec3(0.0f);
		m_Position = glm::vec3(WORLD_SIZE_X / 2, 70, WORLD_SIZE_Z / 2);
		Freefly = false;
		m_isOnGround = false;
		m_AABB.m_Position = m_Position;
	}

	void Player::OnUpdate(GLFWwindow* window, World* world, float dt)
	{
		const float camera_speed = Freefly ? 0.2f : 0.1f;

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		{
			// Take the cross product of the camera's right and up.
			glm::vec3 front = -glm::cross(Camera.GetRight(), Camera.GetUp());
			m_Acceleration += (front * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		{
			glm::vec3 back = glm::cross(Camera.GetRight(), Camera.GetUp());
			m_Acceleration += (back * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		{
			m_Acceleration += (-(Camera.GetRight() * camera_speed));
		}

		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		{
			m_Acceleration += (Camera.GetRight() * camera_speed);
		}

		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		{
			m_Acceleration.y -= camera_speed * 1.35f;
		}

		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		{
			float jump_speed = 0.06f;

			if (!Freefly) 
			{
				if (m_isOnGround)
				{
					m_isOnGround = false;
					m_Acceleration.y += jump_speed * 14.0f;
				}
			}

			else 
			{
				m_Acceleration.y += jump_speed * 3.0f;
			}
		}


		Camera.SetSensitivity(Sensitivity);

		m_Velocity += m_Acceleration;
		m_Acceleration = { 0, 0, 0 };

		// Gravity : 
		if (!Freefly) 
		{
			if (!m_isOnGround) 
			{
				m_Velocity.y -= 0.2 * dt;
			}

			m_isOnGround = false;
		}

		// Test collisions on three axes 
		m_Position.x += m_Velocity.x * dt;
		TestBlockCollision(m_Position, world, glm::vec3(m_Velocity.x, 0.0f, 0.0f));
		m_Position.y += m_Velocity.y * dt;
		TestBlockCollision(m_Position, world, glm::vec3(0.0f, m_Velocity.y, 0.0f));
		m_Position.z += m_Velocity.z * dt;
		TestBlockCollision(m_Position, world, glm::vec3(0.0f, 0.0f, m_Velocity.z));

		// 
		m_AABB.SetPosition(m_Position);

		m_Velocity.x *= 0.825f;
		m_Velocity.z *= 0.825f;

		if (Freefly) {
			m_Velocity.y *= 0.9f;
		}

		Camera.SetPosition(m_Position);
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

	void Player::TestBlockCollision(glm::vec3& position, World* world, glm::vec3 vel)
	{
		for (int x = position.x - m_AABB.m_Dimensions.x; x < position.x + m_AABB.m_Dimensions.x; x++)
		{
			for (int y = position.y - m_AABB.m_Dimensions.y; y < position.y + 0.7; y++)
			{
				for (int z = position.z - m_AABB.m_Dimensions.z; z < position.z + m_AABB.m_Dimensions.z; z++)
				{
					if (x >= 0 && x < WORLD_SIZE_X - 1 &&
						y >= 0 && y < WORLD_SIZE_Y - 1 &&
						z >= 0 && z < WORLD_SIZE_Z - 1)
					{

						Block block = world->GetBlock(x, y, z);

						if (block.block != 0)
						{
							if (vel.y > 0)
							{
								position.y = y - m_AABB.m_Dimensions.y;
								m_Velocity.y = 0;
							}

							else if (vel.y < 0)
							{
								m_isOnGround = true;
								position.y = y + m_AABB.m_Dimensions.y + 1;
								m_Velocity.y = 0;
							}

							if (vel.x > 0)
							{
								position.x = x - m_AABB.m_Dimensions.x;
							}

							else if (vel.x < 0)
							{
								position.x = x + m_AABB.m_Dimensions.x + 1.0f;
							}

							if (vel.z > 0) {
								position.z = z - m_AABB.m_Dimensions.z;
							}
							else if (vel.z < 0) 
							{
								position.z = z + m_AABB.m_Dimensions.z + 1.0f;
							}
						}
					}
				}
			}
		}
	}

	void Player::Jump()
	{
		return;
	}
}