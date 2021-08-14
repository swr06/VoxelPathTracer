#include "BloomRenderer.h"

namespace VoxelRT
{
	namespace BloomRenderer
	{
		static std::unique_ptr<GLClasses::VertexBuffer> BloomFBOVBO;
		static std::unique_ptr<GLClasses::VertexArray> BloomFBOVAO;
		static std::unique_ptr<GLClasses::Framebuffer> BloomAlternateFBO;
		static std::unique_ptr<GLClasses::Shader> BloomMaskShader;
		static std::unique_ptr<GLClasses::Shader> BloomBlurShader;

		void Initialize()
		{
			BloomFBOVBO = std::unique_ptr<GLClasses::VertexBuffer>(new GLClasses::VertexBuffer);
			BloomFBOVAO = std::unique_ptr<GLClasses::VertexArray>(new GLClasses::VertexArray);
			BloomAlternateFBO = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false));
			BloomAlternateFBO->CreateFramebuffer();														 

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

			BloomBlurShader->CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/BloomBlur.glsl");
			BloomBlurShader->CompileShaders();

			BloomMaskShader->CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/BloomMaskFrag.glsl");
			BloomMaskShader->CompileShaders();
		}

		void BlurBloomMip(BloomFBO& bloomfbo, int mip_num, GLuint source_tex, GLuint bright_tex, bool hq)
		{
			GLClasses::Shader& GaussianBlur = *BloomBlurShader;

			GaussianBlur.Use();

			bloomfbo.BindMip(mip_num);

			GaussianBlur.SetInteger("u_Texture", 0);
			GaussianBlur.SetInteger("u_Lod", mip_num);
			GaussianBlur.SetBool("u_HQ", hq);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, bright_tex);

			BloomFBOVAO->Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			BloomFBOVAO->Unbind();
		}

		void RenderBloom(BloomFBO& bloom_fbo, GLuint source_tex, GLuint bright_tex, bool hq)
		{
			// Render the bright parts to a texture
			GLClasses::Shader& BloomBrightShader = *BloomMaskShader;

			BloomAlternateFBO->SetSize(bloom_fbo.GetWidth() * floor(bloom_fbo.m_MipScales[0] * 2.0f), bloom_fbo.GetHeight() * floor(bloom_fbo.m_MipScales[0] * 2.0f));

			BloomBrightShader.Use();
			BloomAlternateFBO->Bind();
			BloomBrightShader.SetInteger("u_Texture", 0);
			BloomBrightShader.SetInteger("u_EmissiveTexture", 1);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, source_tex);

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, bright_tex);

			BloomFBOVAO->Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			BloomFBOVAO->Unbind();

			glDisable(GL_DEPTH_TEST);
			glDisable(GL_CULL_FACE);

			// mip gen

			bool genmip = false;

			if (genmip) {
				glBindTexture(GL_TEXTURE_2D, BloomAlternateFBO->GetTexture(0));
				glGenerateMipmap(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, 0);
			}

			// Blur the mips					
			BlurBloomMip(bloom_fbo, 0, source_tex, BloomAlternateFBO->GetTexture(0), hq);
			BlurBloomMip(bloom_fbo, 1, source_tex, bloom_fbo.m_Mips[0], hq);
			BlurBloomMip(bloom_fbo, 2, source_tex, bloom_fbo.m_Mips[1], hq);
			BlurBloomMip(bloom_fbo, 3, source_tex, bloom_fbo.m_Mips[2], hq);

			return;
		}

		void RecompileShaders()
		{
			BloomBlurShader->Recompile();
			BloomMaskShader->Recompile();
		}
	}
}
