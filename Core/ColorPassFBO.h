#pragma once

#include <glad/glad.h>

#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "../Core/Application/Logger.h"

namespace VoxelRT
{
	class ColorPassFBO
	{
	public:

		ColorPassFBO(uint32_t w = 32, uint32_t h = 32);
		~ColorPassFBO();

		ColorPassFBO(const ColorPassFBO&) = delete;
		ColorPassFBO operator=(ColorPassFBO const&) = delete;

		ColorPassFBO& operator=(ColorPassFBO&& other)
		{
			std::swap(*this, other);
			return *this;
		}

		ColorPassFBO(ColorPassFBO&& v)
		{
			m_Width = v.m_Width;
			m_Height = v.m_Height;

			m_ColorTexture = v.m_ColorTexture;
			m_NormalTexture = v.m_NormalTexture;
			m_PBR = v.m_PBR;
			m_FBO = v.m_FBO;

			v.m_ColorTexture = 0;
			v.m_NormalTexture = 0;
			v.m_PBR = 0;
			v.m_FBO = 0;
		}

		void Bind() const
		{
			glBindFramebuffer(GL_FRAMEBUFFER, m_FBO);
			glViewport(0, 0, m_Width, m_Height);
		}

		void Unbind() const
		{
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}

		void SetDimensions(uint32_t w, uint32_t h);

		std::pair<uint32_t, uint32_t> GetDimensions()
		{
			return { m_Width, m_Height };
		}

		uint32_t GetWidth() { return m_Width; }
		uint32_t GetHeight() { return m_Height; }

		inline GLuint GetFramebufferID() const noexcept { return m_FBO; }
		inline GLuint GetColorTexture() const noexcept { return m_ColorTexture; }
		inline GLuint GetNormalTexture() const noexcept { return m_NormalTexture; }
		inline GLuint GetPBRTexture() const noexcept { return m_PBR; }
		void GenerateFramebuffers();

	private:

		void DeleteEverything();

		GLuint m_ColorTexture = 0;
		GLuint m_NormalTexture = 0;
		GLuint m_PBR = 0;

		GLuint m_FBO = 0;
		uint32_t m_Width = 0;
		uint32_t m_Height = 0;
	};
}