#include "Texture3D.h"

VoxelRT::Texture3D::Texture3D()
{
	m_ID = 0;
}

void VoxelRT::Texture3D::CreateTexture(int w, int h, int d, void* data)
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
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, w, h, d, 0, GL_RED, GL_UNSIGNED_BYTE, data);
}
