#pragma once

#include <glad/glad.h>

#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "../Application/Logger.h"

namespace Clouds
{
	class CloudFBO
	{
	public:

		CloudFBO(uint32_t w = 32, uint32_t h = 32);
		~CloudFBO();

		CloudFBO(const CloudFBO&) = delete;
		CloudFBO operator=(CloudFBO const&) = delete;

		CloudFBO& operator=(CloudFBO&& other)
		{
			std::swap(*this, other);
			return *this;
		}

		CloudFBO(CloudFBO&& v)
		{
			m_Width = v.m_Width;
			m_Height = v.m_Height;

			m_CloudData = v.m_CloudData;
			m_FBO = v.m_FBO;

			v.m_CloudData = 0;
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
		inline GLuint GetCloudTexture() const noexcept { return m_CloudData; }
		void GenerateFramebuffers();

	private:

		void DeleteEverything();

		GLuint m_CloudData = 0;

		GLuint m_FBO = 0;
		uint32_t m_Width = 0;
		uint32_t m_Height = 0;
	};
}