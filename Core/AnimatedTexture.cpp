#include "AnimatedTexture.h"
#include <string>
#include <sstream>

#include "../Core/GLClasses/stb_image.h"

void VoxelRT::AnimatedTexture::Create(const std::string& path, int res, int frames)
{
	glGenTextures(1, &m_ID);
	glBindTexture(GL_TEXTURE_3D, m_ID);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, res, res, frames+1, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);

	for (int i = 0; i <= frames; i++) {
		int w, h, bpp;
		std::string current_path = path;
		current_path += std::string("/") + std::to_string(i) + std::string(".png");
		unsigned char* image = stbi_load(current_path.c_str(), &w, &h, &bpp, 0);
		glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, i, res, res, 1, GL_RGB, GL_UNSIGNED_BYTE, image);
	}

	glBindTexture(GL_TEXTURE_3D, 0);
}
