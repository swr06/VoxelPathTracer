#pragma once

#include <glad/glad.h>

#include <iostream>
#include <string>

namespace Clouds
{
	class NoiseTexture3D
	{
	public :
		NoiseTexture3D();
		void CreateTexture(int w, int h, int d, void* data);
		void Delete();

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
