#pragma once

#include <glad/glad.h>

#include <array>
#include <vector>
#include <utility>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLClasses/VertexBuffer.h"
#include "GLClasses/VertexArray.h"
#include "GLClasses/IndexBuffer.h"
#include "GLClasses/Shader.h"
#include "GLClasses/Texture.h"

#include "OrthographicCamera.h"

namespace VoxelRT
{
	class Renderer2D
	{
	public:

		Renderer2D();
	
		void RenderQuad(const glm::vec2& position, GLClasses::Texture* texture, OrthographicCamera* camera);
		void RenderQuad(const glm::vec2& position, GLClasses::Texture* texture, OrthographicCamera* camera, int w, int h);

	private : 

		GLClasses::VertexBuffer m_VBO;
		GLClasses::VertexArray m_VAO;
		GLClasses::IndexBuffer m_IBO;
		GLClasses::Shader m_DefaultShader;
	};
}