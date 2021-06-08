#include "NoiseRenderer.h"
#include <thread>

namespace Clouds
{
	void RenderNoise(NoiseTexture3D& tex, int slices, bool detail)
	{
		GLuint FBO = 0;
		GLClasses::Shader NoiseShader;

		float Vertices[] =
		{
			-1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
			 1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
			 1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
		};

		GLClasses::VertexBuffer VBO;
		GLClasses::VertexArray VAO;
		VAO.Bind();
		VBO.Bind();
		VBO.BufferData(sizeof(Vertices), Vertices, GL_STATIC_DRAW);
		VBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
		VBO.VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
		VAO.Unbind();

		if (detail)
		{
			NoiseShader.CreateShaderProgramFromFile("Core/Shaders/Clouds/FBOVert.glsl", "Core/Shaders/Clouds/NoiseDetailFrag.glsl");
		}

		else
		{
			NoiseShader.CreateShaderProgramFromFile("Core/Shaders/Clouds/FBOVert.glsl", "Core/Shaders/Clouds/NoiseFrag.glsl");
		}

		NoiseShader.CompileShaders();

		glGenFramebuffers(1, &FBO);
		glBindFramebuffer(GL_FRAMEBUFFER, FBO);

		for (int i = 0; i < slices; i++)
		{
			glViewport(0, 0, tex.GetWidth(), tex.GetHeight());
			glBindFramebuffer(GL_FRAMEBUFFER, FBO);
			glFramebufferTexture3D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_3D, tex.GetTextureID(), 0, i);
			glClear(GL_COLOR_BUFFER_BIT);

			NoiseShader.Use();
			NoiseShader.SetFloat("u_CurrentSlice", (float)i / (float)slices);
			NoiseShader.SetVector2f("u_Dims", glm::vec2(tex.GetWidth(), tex.GetHeight()));
			
			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			std::this_thread::sleep_for(std::chrono::milliseconds(4));;
		}
	}
}