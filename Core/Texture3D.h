#pragma once

#include <glad/glad.h>

#include <iostream>
#include <string>

namespace VoxelRT
{
	class Texture3D
	{
	public :
		Texture3D();
		void CreateTexture(int w, int h, int d, void* data);

		inline int GetWidth() const { return m_Width; }
		inline int GetHeight() const { return m_Height; }
		inline int GetDepth() const { return m_Depth; }
		inline GLuint GetTextureID() const { return m_ID; }

	private :
		GLuint m_ID;
		int m_Width;
		int m_Height;
		int m_Depth;
	};
}
