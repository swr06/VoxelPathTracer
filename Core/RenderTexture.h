#pragma once

#include <iostream>
#include <glad/glad.h>
#include <string>
#include <vector>

namespace GLClasses
{
	class RenderTexture
	{
	public:
		RenderTexture() : m_ID(0), m_Width(1), m_Height(1) {}
		~RenderTexture();

		RenderTexture(const RenderTexture&) = delete;
		RenderTexture operator=(RenderTexture const&) = delete;
		RenderTexture(RenderTexture&& v)
		{
			m_ID = v.m_ID;
			m_Width = v.m_Width;
			m_Height = v.m_Height;

			v.m_ID = 0;
		}

		void Resize(int w, int h)
		{
			if (m_Width != w || m_Height != h)
			{
				glDeleteTextures(1, &m_ID);
				m_ID = 0;

				CreateRenderTexture(w, h);
			}
		}

		void CreateRenderTexture(int, int);

		GLuint GetID() const noexcept
		{
			return m_ID;
		}

		int GetWidth() const noexcept { return m_Width; }
		int GetHeight() const noexcept { return m_Height; }

	private:

		GLuint m_ID = 0;
		int m_Width = 0;
		int m_Height = 0;
	};
}