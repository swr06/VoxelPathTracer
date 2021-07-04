#include "Framebuffer.h"

#ifndef NULL
	#define NULL 0
#endif

namespace GLClasses
{
    Framebuffer::Framebuffer(unsigned int w, unsigned int h, std::vector<FORMAT> format, bool has_depth_attachment) :
        m_FBO(0), m_FBWidth(w), m_FBHeight(h), m_HasDepthMap(has_depth_attachment), m_Format(format)
    {
       // CreateFramebuffer();
    }

    Framebuffer::Framebuffer(unsigned int w, unsigned int h, FORMAT format, bool has_depth_attachment) :
        m_FBO(0), m_FBWidth(w), m_FBHeight(h), m_HasDepthMap(has_depth_attachment), m_Format({ format })
    {
        // CreateFramebuffer();
    }

	Framebuffer::~Framebuffer()
	{
        glDeleteFramebuffers(1, &m_FBO);

        for (int i = 0; i < m_TextureAttachments.size(); i++)
        {
            glDeleteTextures(1, &m_TextureAttachments[i]);
            m_TextureAttachments[i];
        }
        glDeleteRenderbuffers(1, &m_DepthStencilBuffer);
	}

	void Framebuffer::CreateFramebuffer()
	{
        std::vector<GLenum> DrawBuffers;

        if (m_Format.size() > 31)
        {
            throw "Too many attachments.";
        }

        int w = m_FBWidth;
        int h = m_FBHeight;

        glGenFramebuffers(1, &m_FBO);
        glBindFramebuffer(GL_FRAMEBUFFER, m_FBO);


        for (int i = 0; i < m_Format.size(); i++)
        {
            GLenum format = m_Format[i].Format;
            GLenum type = m_Format[i].Type;
            bool min = m_Format[i].MinFilter;
            bool mag = m_Format[i].MagFilter;
            m_TextureAttachments.emplace_back();
            glGenTextures(1, &m_TextureAttachments.at(i));
            glBindTexture(GL_TEXTURE_2D, m_TextureAttachments.at(i));
            glTexImage2D(GL_TEXTURE_2D, 0, format, w, h, 0, m_Format[i].InternalFormat, type, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min ? GL_LINEAR : GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag ? GL_LINEAR : GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + i, GL_TEXTURE_2D, m_TextureAttachments.at(i), 0);
            DrawBuffers.push_back(GL_COLOR_ATTACHMENT0+i);
        }

        glDrawBuffers(m_Format.size(), DrawBuffers.data());
       
        if (m_HasDepthMap)
        {
            glGenRenderbuffers(1, &m_DepthStencilBuffer);
            glBindRenderbuffer(GL_RENDERBUFFER, m_DepthStencilBuffer);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, w, h);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, m_DepthStencilBuffer);
        }
       
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            VoxelRT::Logger::Log("Fatal error! Framebuffer creation failed!");
        }

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

    
}