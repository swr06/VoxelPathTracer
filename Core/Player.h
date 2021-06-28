#pragma once

#include <glm/glm.hpp>
#include "World.h"
#include "FpsCamera.h"

#include <glfw/glfw3.h>

namespace VoxelRT
{
	class Player
	{
	public :

		Player();
		bool OnUpdate(GLFWwindow* window, World* world, float dt);

		FPSCamera Camera;
		bool TestBlockCollision(const glm::vec3& position, World* world);

		bool InWater = false;
		bool Freefly = false;
		float Sensitivity = 0.25;
		float Speed = 0.045f;
	};
}