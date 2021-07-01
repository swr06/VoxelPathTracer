#pragma once

#include <iostream>
#include <string>
#include <glad/glad.h>
#include "../Application/Logger.h"

namespace GLClasses
{
	class Framebuffer
	{
	public :
		Framebuffer(unsigned int w = 16, unsigned int h = 16, bool hdr = true, bool has_depth_attachment = false, bool linear = true);
		~Framebuffer();

		Framebuffer(const Framebuffer&) = delete;
		Framebuffer operator=(Framebuffer const&) = delete;

		Framebuffer& operator=(Framebuffer&& other)
		{
			std::swap(*this, other);
			return *this;
		}

		Framebuffer(Framebuffer&& v) : m_IsHDR(v.m_IsHDR), m_HasDepthMap(v.m_HasDepthMap)
		{
			m_FBO = v.m_FBO;
			m_TextureAttachment = v.m_TextureAttachment;
			m_DepthStencilBuffer = v.m_DepthStencilBuffer;
			m_FBWidth = v.m_FBWidth;
			m_FBHeight = v.m_FBHeight;

			v.m_FBO = 0;
			v.m_TextureAttachment = 0;
			v.m_DepthStencilBuffer = 0;
		}

		void Bind() const
		{
			glBindFramebuffer(GL_FRAMEBUFFER, m_FBO);
			glViewport(0, 0, m_FBWidth, m_FBHeight);
		}


		void Unbind() const
		{
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}

		void SetSize(uint32_t width, uint32_t height)
		{
			if (width != m_FBWidth || height != m_FBHeight)
			{
				glDeleteFramebuffers(1, &m_FBO);
				glDeleteTextures(1, &m_TextureAttachment);
				glDeleteRenderbuffers(1, &m_DepthStencilBuffer);

				m_FBO = 0;
				m_TextureAttachment = 0;
				m_DepthStencilBuffer = 0;
				m_FBWidth = width;
				m_FBHeight = height;
				CreateFramebuffer();
			}
		}

		inline GLuint GetTexture() const 
		{
			return m_TextureAttachment;
		}

		inline GLuint GetDepthStencilBuffer() const
		{
			return m_DepthStencilBuffer;
		}

		inline GLuint GetFramebuffer() const noexcept { return m_FBO; }
		inline unsigned int GetWidth() const noexcept { return m_FBWidth; }
		inline unsigned int GetHeight() const noexcept { return m_FBHeight; }

		inline bool IsHDR() const
		{
			return m_IsHDR;
		}

		void SetExposure(float exp)
		{
			if (!m_IsHDR)
			{
				throw "SetExposure(float) called on a Framebuffer that is not HDR!";
			}

			m_Exposure = exp;
		}

		float GetExposure() const noexcept
		{
			return m_Exposure;
		}

		// Creates the framebuffer with the appropriate settings
		void CreateFramebuffer();

	private :

		GLuint m_FBO; // The Framebuffer object
		GLuint m_TextureAttachment; // The actual texture attachment
		GLuint m_DepthStencilBuffer;
		int m_FBWidth;
		int m_FBHeight;
		const bool m_IsHDR;
		const bool m_HasDepthMap;
		float m_Exposure = 0.0f;
		bool m_Linear;
	};
}