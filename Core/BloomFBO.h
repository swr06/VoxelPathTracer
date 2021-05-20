#pragma once

#include <iostream>
#include <glad/glad.h>

namespace VoxelRT
{
	class BloomFBO
	{
	public :
		BloomFBO(int w, int h);
		~BloomFBO();

		BloomFBO(const BloomFBO&) = delete;
		BloomFBO operator=(BloomFBO const&) = delete;

		BloomFBO(BloomFBO&& v)
		{
			m_Framebuffer = v.m_Framebuffer;
			m_Mip0 = v.m_Mip0;
			m_Mip1 = v.m_Mip1;
			m_w = v.m_w;
			m_h = v.m_h;

			v.m_Framebuffer = 0;
			v.m_Mip0 = 0;
			v.m_Mip1 = 0;
			v.m_w = -1;
			v.m_h = -1;
		}

		void SetSize(int w, int h)
		{
			if (w != m_w || h != m_h)
			{
				DeleteEverything();
				Create(w, h);
				m_w = w;
				m_h = h;
			}
		}

		GLuint m_Framebuffer;
		GLuint m_Mip0;
		GLuint m_Mip1;
		GLuint m_Mip2;
		GLuint m_Mip3;

		inline int GetWidth() const { return m_w; }
		inline int GetHeight() const { return m_h; }

		const float m_mipscale1 = 0.9f;
		const float m_mipscale2 = 0.5f;
		const float m_mipscale3 = 0.3f;
		const float m_mipscale4 = 0.15f;

	private :
		void DeleteEverything();
		void Create(int w, int h);

		int m_w = -1, m_h = -1;
	};
}