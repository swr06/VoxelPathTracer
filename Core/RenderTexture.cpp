#include "RenderTexture.h"

GLClasses::RenderTexture::~RenderTexture()
{
	glDeleteTextures(1, &m_ID);
	m_ID = 0;
}

void GLClasses::RenderTexture::CreateRenderTexture(int w, int h)
{
	if (m_ID != 0)
	{
		throw "Trying to create a texture when one is already created!";
	}

	glGenTextures(1, &m_ID);
	glBindTexture(GL_TEXTURE_2D, m_ID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, NULL);

	glBindTexture(GL_TEXTURE_2D, 0);

	m_Width = w;
	m_Height = h;
}
