#pragma once

#include <glad/glad.h>

#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "../Core/Application/Logger.h"

namespace VoxelRT
{
	class InitialRTFBO
	{
	public:

		InitialRTFBO(uint32_t w = 32, uint32_t h = 32);
		~InitialRTFBO();

		InitialRTFBO(const InitialRTFBO&) = delete;
		InitialRTFBO operator=(InitialRTFBO const&) = delete;

		InitialRTFBO& operator=(InitialRTFBO&& other)
		{
			std::swap(*this, other);
			return *this;
		}

		InitialRTFBO(InitialRTFBO&& v)
		{
			m_Width = v.m_Width;
			m_Height = v.m_Height;

			m_PositionTexture = v.m_PositionTexture;
			m_NormalTexture = v.m_NormalTexture;
			m_DataTexture = v.m_DataTexture;
			m_FBO = v.m_FBO;

			v.m_PositionTexture = 0;
			v.m_NormalTexture = 0;
			v.m_DataTexture = 0;
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
		inline GLuint GetPositionTexture() const noexcept { return m_PositionTexture; }
		inline GLuint GetNormalTexture() const noexcept { return m_NormalTexture; }
		inline GLuint GetDataTexture() const noexcept { return m_DataTexture; }
		void GenerateFramebuffers();

	private:

		void DeleteEverything();

		GLuint m_PositionTexture = 0;
		GLuint m_NormalTexture = 0;
		GLuint m_DataTexture = 0;

		GLuint m_FBO = 0;
		uint32_t m_Width = 0;
		uint32_t m_Height = 0;
	};
}