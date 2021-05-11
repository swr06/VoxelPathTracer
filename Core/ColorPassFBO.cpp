#include "ColorPassFBO.h"

namespace VoxelRT
{
	ColorPassFBO::ColorPassFBO(uint32_t w, uint32_t h)
	{
		m_Width = w;
		m_Height = h;
	}

	ColorPassFBO::~ColorPassFBO()
	{
		DeleteEverything();
	}

	void ColorPassFBO::SetDimensions(uint32_t w, uint32_t h)
	{
		if (w != m_Width || h != m_Height)
		{
			DeleteEverything();
			m_Width = w;
			m_Height = h;
			GenerateFramebuffers();
		}
	}

	void ColorPassFBO::GenerateFramebuffers()
	{
		GLenum DrawBuffers[3] = {
			GL_COLOR_ATTACHMENT0,
			GL_COLOR_ATTACHMENT1,
			GL_COLOR_ATTACHMENT2
		};

		glGenFramebuffers(1, &m_FBO);
		glBindFramebuffer(GL_FRAMEBUFFER, m_FBO);

		// Create an HDR color buffer
		glGenTextures(1, &m_ColorTexture);
		glBindTexture(GL_TEXTURE_2D, m_ColorTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, m_Width, m_Height, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_ColorTexture, 0);

		// Normal texture
		glGenTextures(1, &m_NormalTexture);
		glBindTexture(GL_TEXTURE_2D, m_NormalTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, m_Width, m_Height, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, m_NormalTexture, 0);

		// data texture
		glGenTextures(1, &m_PBR);
		glBindTexture(GL_TEXTURE_2D, m_PBR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, m_Width, m_Height, 0, GL_RGBA, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, m_PBR, 0);

		glDrawBuffers(3, DrawBuffers);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		{
			VoxelRT::Logger::Log("Geometry Render buffer is not complete!\n\tWIDTH : " +
				std::to_string(m_Width) + "\n\tHEIGHT : " + std::to_string(m_Height));
		}

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	void ColorPassFBO::DeleteEverything()
	{
		glDeleteFramebuffers(1, &m_FBO);
		glDeleteTextures(1, &m_ColorTexture);
		glDeleteTextures(1, &m_NormalTexture);
		glDeleteTextures(1, &m_PBR);
		m_FBO = 0;
		m_ColorTexture = 0;
		m_NormalTexture = 0;
		m_PBR = 0;
	}
}