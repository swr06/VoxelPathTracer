#include "BloomRenderer.h"

namespace VoxelRT
{
	namespace BloomRenderer
	{
		static std::unique_ptr<GLClasses::VertexBuffer> BloomFBOVBO;
		static std::unique_ptr<GLClasses::VertexArray> BloomFBOVAO;
		static std::unique_ptr<GLClasses::Framebuffer> BloomAlternateFBO[4];
		static std::unique_ptr<GLClasses::Shader> BloomMaskShader;
		static std::unique_ptr<GLClasses::Shader> BloomBlurShader;

		void Initialize()
		{
			BloomFBOVBO = std::unique_ptr<GLClasses::VertexBuffer>(new GLClasses::VertexBuffer);
			BloomFBOVAO = std::unique_ptr<GLClasses::VertexArray>(new GLClasses::VertexArray);
			BloomAlternateFBO[0] = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO[0]->CreateFramebuffer();
			BloomAlternateFBO[1] = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO[1]->CreateFramebuffer();
			BloomAlternateFBO[2] = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO[2]->CreateFramebuffer();
			BloomAlternateFBO[3] = std::unique_ptr<GLClasses::Framebuffer>(new GLClasses::Framebuffer(64, 64, true));
			BloomAlternateFBO[3]->CreateFramebuffer();

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
			int w = floor(bloomfbo.GetWidth() * bloomfbo.m_MipScales[mip_num]), h = floor(bloomfbo.GetHeight() * bloomfbo.m_MipScales[mip_num]);
			std::unique_ptr<GLClasses::Framebuffer>* fbo = &BloomAlternateFBO[mip_num];

			GLClasses::Shader& GaussianBlur = *BloomBlurShader;
			GLClasses::Shader& BloomBrightShader = *BloomMaskShader;


			fbo->get()->SetSize(w, h);

			BloomBrightShader.Use();
			fbo->get()->Bind();

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

			bloomfbo.BindMip(mip_num);

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
