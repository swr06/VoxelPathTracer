#pragma once

#include <glm/glm.hpp>
#include "World.h"
#include "FpsCamera.h"

#include <glfw/glfw3.h>

#include "Physics/AABB.h"
#include "SoundManager.h"

namespace VoxelRT
{
	class Player
	{
	public :

		Player();
		void OnUpdate(GLFWwindow* window, World* world, float dt, int frame);

		FPSCamera Camera;
		void TestBlockCollision(glm::vec3& position, World* world, glm::vec3 vel);
		void Jump();

		bool InWater = false;
		bool Freefly = false;
		float Sensitivity = 0.25;
		float Speed = 0.045f;


		glm::vec3 m_Position;
		glm::vec3 m_Velocity;
		glm::vec3 m_Acceleration;
		AABB m_AABB;
		bool m_isOnGround;
		bool DisableCollisions = false;

	private :


	protected :



	};
}