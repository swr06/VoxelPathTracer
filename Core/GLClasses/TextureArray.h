#pragma once

#include <iostream>
#include <string>
#include <array>
#include <unordered_map>
#include <set>
#include <algorithm>

#include "stb_image.h"
#include <glad/glad.h>
#include <glfw/glfw3.h>

namespace GLClasses
{
	class TextureArray
	{
	public :

		TextureArray();
		~TextureArray()
		{
			glDeleteTextures(1, &m_TextureArray);
		}

		TextureArray(const TextureArray&) = delete;
		TextureArray operator=(TextureArray const&) = delete;
		TextureArray(TextureArray&& v)
		{
			m_TextureLocations = v.m_TextureLocations;
			m_TextureArray = v.m_TextureArray;
			v.m_TextureArray = 0;
		}

		void CreateArray(std::vector<std::string> paths, std::pair<int, int> texture_size, bool is_srgb ,bool gen_mips = true, GLint mag_filter = GL_NEAREST, bool limit_textures = false);
		GLuint GetTextureArray() const noexcept { return m_TextureArray; }

		int GetTexture(const std::string& s) const 
		{ 
			return m_TextureLocations.at(s); // Todo : Make this secure
		}

	private : 

		GLuint m_TextureArray = 0;
		std::unordered_map<std::string, int> m_TextureLocations;
	};
}
