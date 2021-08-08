#include "NoiseTexture3D.h"

Clouds::NoiseTexture3D::NoiseTexture3D()
{
	m_ID = 0;
}

void Clouds::NoiseTexture3D::CreateTexture(int w, int h, int d, void* data)
{
	if (m_ID > 0)
	{
		glDeleteTextures(1, &m_ID);
		m_ID = 0;
	}

	m_Width = w;
	m_Height = h;
	m_Depth = d;

	glGenTextures(1, &m_ID);
	glBindTexture(GL_TEXTURE_3D, m_ID);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);

	//glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA16F, w, h, d, 0, GL_RGBA, GL_FLOAT, data);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, w, h, d, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
}

void Clouds::NoiseTexture3D::Delete()
{
	glDeleteTextures(1, &m_ID);
}
