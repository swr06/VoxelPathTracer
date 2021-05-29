#include "BloomRenderer.h"

namespace VoxelRT
{
	namespace BloomRenderer
	{
		static std::unique_ptr<GLClasses::VertexBuffer> BloomFBOVBO;
		static std::unique_ptr<GLClasses::VertexArray> BloomFBOVAO;
		static std::unique_ptr<GLClasses::Framebuffer> BloomAlternateFBO;
		static std::unique_ptr<GLClasses::Framebuffer> BloomAlternateFBO1;
		static std::unique_ptr<GLClasses::Framebuffer> BloomAlternateFBO2;
		static std::unique_ptr<GLClasses::Framebuffer> BloomAlternateFBO3;
		static std::unique_ptr<GLClasses::Shader> BloomMaskShader;
		static std::unique_ptr<GLClasses::Shader> BloomBlurShader;

		void Initialize()
		{
			BloomFBOVBO = std::unique_ptr<GLClasses::VertexBuffer>(new GLClasses::VertexBuffer);
			BloomFBOVAO = std::unique_ptr<GLClasses::VertexArray>(new GLClasses::VertexArray);
			BloomAlternateFBO = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO->CreateFramebuffer();
			BloomAlternateFBO1 = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO1->CreateFramebuffer();
			BloomAlternateFBO2 = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO2->CreateFramebuffer();
			BloomAlternateFBO3 = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO3->CreateFramebuffer();
			BloomBlurShader = std::unique_ptr<GLClasses::Shader>(new GLClasses::Shader);
			BloomMaskShader = std::unique_ptr<GLClasses::Shader>(new GLClasses::Shader);

			float QuadVertices[] =
			{
				-1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
				 1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
				 1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
			};

			BloomFBOVBO->BufferData(sizeof(QuadVertices), QuadVertices, GL_STATIC_DRAW);
			BloomFBOVAO->Bind();
			BloomFBOVBO->Bind();
			BloomFBOVBO->VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
			BloomFBOVBO->VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
			BloomFBOVAO->Unbind();

			BloomBlurShader->CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/Gaussian5TapSinglePass.glsl");
			BloomBlurShader->CompileShaders();

			BloomMaskShader->CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/BloomMaskFrag.glsl");
			BloomMaskShader->CompileShaders();
		}

		void BlurBloomMip(BloomFBO& bloomfbo, int mip_num, GLuint source_tex, GLuint bright_tex)
		{
			GLenum buffer;
			int w, h;
			float mip = 0.0f;
			std::unique_ptr<GLClasses::Framebuffer>* fbo;

			GLClasses::Shader& GaussianBlur = *BloomBlurShader;
			GLClasses::Shader& BloomBrightShader = *BloomMaskShader;

			switch (mip_num)
			{
			case 0:
			{
				buffer = GL_COLOR_ATTACHMENT0;
				w = floor(bloomfbo.m_mipscale1 * bloomfbo.GetWidth());
				h = floor(bloomfbo.m_mipscale1 * bloomfbo.GetHeight());
				mip = 0.0f;
				fbo = &BloomAlternateFBO;

				break;
			}

			case 1:
			{
				buffer = GL_COLOR_ATTACHMENT1;
				w = floor(bloomfbo.m_mipscale2 * bloomfbo.GetWidth());
				h = floor(bloomfbo.m_mipscale2 * bloomfbo.GetHeight());
				mip = 0.0f;
				fbo = &BloomAlternateFBO1;

				break;
			}

			case 2:
			{
				buffer = GL_COLOR_ATTACHMENT2;
				w = floor(bloomfbo.m_mipscale3 * bloomfbo.GetWidth());
				h = floor(bloomfbo.m_mipscale3 * bloomfbo.GetHeight());
				mip = 0.0f;
				fbo = &BloomAlternateFBO2;

				break;
			}

			case 3:
			{
				buffer = GL_COLOR_ATTACHMENT3;
				w = floor(bloomfbo.m_mipscale4 * bloomfbo.GetWidth());
				h = floor(bloomfbo.m_mipscale4 * bloomfbo.GetHeight());
				mip = 0.0f;
				fbo = &BloomAlternateFBO3;

				break;
			}

			default:
			{
				throw "INVALID VALUE PASSED TO BLOOM RENDERER!";
			}
			}

			fbo->get()->SetSize(w, h);

			BloomBrightShader.Use();
			fbo->get()->Bind();
			glViewport(0, 0, w, h);

			BloomBrightShader.SetInteger("u_Texture", 0);
			BloomBrightShader.SetInteger("u_EmissiveTexture", 1);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, source_tex);

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, bright_tex);

			BloomFBOVAO->Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			BloomFBOVAO->Unbind();

			///////////////

			GaussianBlur.Use();

			glBindFramebuffer(GL_FRAMEBUFFER, bloomfbo.m_Framebuffer);
			glDrawBuffer(buffer);
			glViewport(0, 0, w, h);

			GaussianBlur.SetInteger("u_Texture", 0);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, fbo->get()->GetTexture());

			BloomFBOVAO->Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			BloomFBOVAO->Unbind();
		}

		void RenderBloom(BloomFBO& bloom_fbo, GLuint source_tex, GLuint bright_tex)
		{
			// Render the bright parts to a texture

			glDisable(GL_DEPTH_TEST);
			glDisable(GL_CULL_FACE);

			// Blur the mips					
			BlurBloomMip(bloom_fbo, 0, source_tex, bright_tex);
			BlurBloomMip(bloom_fbo, 1, source_tex, bright_tex);
			BlurBloomMip(bloom_fbo, 2, source_tex, bright_tex);
			BlurBloomMip(bloom_fbo, 3, source_tex, bright_tex);

			return;
		}

		void RecompileShaders()
		{
			BloomBlurShader->Recompile();
			BloomMaskShader->Recompile();
		}
	}
}
