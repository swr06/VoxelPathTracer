#pragma once

#include <iostream>
#include <string>
#include <array>

#include "GLClasses/CubeTextureMap.h"
#include "GLClasses/Shader.h"
#include "GLClasses/VertexBuffer.h"
#include "GLClasses/VertexArray.h"
#include "FpsCamera.h"
#include "AtmosphereRenderCubemap.h"

namespace VoxelRT
{
	class AtmosphereRenderer
	{
	public:
		AtmosphereRenderer();
		void RenderAtmosphere(AtmosphereRenderMap& map, const glm::vec3& sun_direction, int steps, int lsteps);
		void Recompile();

	private :
		GLClasses::VertexBuffer m_VBO;
		GLClasses::VertexArray m_VAO;
		GLClasses::Shader m_AtmosphereShader;
	};
}