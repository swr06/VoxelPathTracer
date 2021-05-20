#include "BloomFBO.h"

namespace VoxelRT
{ 
	BloomFBO::BloomFBO(int w, int h)
	{
		Create(w, h);
	}

	void BloomFBO::Create(int w, int h)
	{
		m_w = w;
		m_h = h;

		int w0, h0, w1, h1, w2, h2, w3, h3;

		w0 = w * m_mipscale1;
		h0 = h * m_mipscale1;

		w1 = w * m_mipscale2;
		h1 = h * m_mipscale2;

		w2 = w * m_mipscale3;
		h2 = h * m_mipscale3;

		w3 = w * m_mipscale4;
		h3 = h * m_mipscale4;

		glGenTextures(1, &m_Mip0);
		glBindTexture(GL_TEXTURE_2D, m_Mip0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w0, h0, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

		glGenTextures(1, &m_Mip1);
		glBindTexture(GL_TEXTURE_2D, m_Mip1);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w1, h1, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

		glGenTextures(1, &m_Mip2);
		glBindTexture(GL_TEXTURE_2D, m_Mip2);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w2, h2, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

		glGenTextures(1, &m_Mip3);
		glBindTexture(GL_TEXTURE_2D, m_Mip3);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w3, h3, 0, GL_RGB, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);


		glGenFramebuffers(1, &m_Framebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, m_Framebuffer);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_Mip0, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, m_Mip1, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, m_Mip2, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, m_Mip3, 0);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		{
			std::cout << "\nBloom fbo wasn't created\n";
		}

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	BloomFBO::~BloomFBO()
	{
		DeleteEverything();
	}

	void BloomFBO::DeleteEverything()
	{
		glDeleteTextures(1, &m_Mip0);
		glDeleteTextures(1, &m_Mip1);
		glDeleteFramebuffers(1, &m_Framebuffer);

		m_Mip0 = 0;
		m_Mip1 = 0;
		m_Framebuffer = 0;
	}
}