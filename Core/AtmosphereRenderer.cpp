#include "AtmosphereRenderer.h"

#include <glm/glm.hpp>

namespace VoxelRT
{
    AtmosphereRenderer::AtmosphereRenderer() : m_VBO(GL_ARRAY_BUFFER)
	{
        m_AtmosphereShader.CreateShaderProgramFromFile("Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/AtmosphereFrag.glsl");
        m_AtmosphereShader.CompileShaders();

        float QuadVertices[] =
        {
            -1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
             1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
             1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
        };

        m_VBO.BufferData(sizeof(QuadVertices), QuadVertices, GL_STATIC_DRAW);
        m_VAO.Bind();
        m_VBO.Bind();
        m_VBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
        m_VBO.VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
        m_VAO.Unbind();
	}

    void AtmosphereRenderer::RenderAtmosphere(AtmosphereRenderMap& map, const glm::vec3& sun_direction, int steps, int lsteps)
    {
        glDepthMask(GL_FALSE);
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);

        glm::mat4 projection_matrix = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, 600.0f);
        glm::vec3 center = glm::vec3(0.0f, 5000.0f, 0.0f);

        std::array<glm::mat4, 6> view_matrices =
        {
            glm::lookAt(center, center + glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
            glm::lookAt(center, center + glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
            glm::lookAt(center, center + glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f,  0.0f,  1.0f)),
            glm::lookAt(center, center + glm::vec3(0.0f,-1.0f, 0.0f), glm::vec3(0.0f,  0.0f, -1.0f)),
            glm::lookAt(center, center + glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
            glm::lookAt(center, center + glm::vec3(0.0f, 0.0f,-1.0f), glm::vec3(0.0f, -1.0f,  0.0f))
        };

        m_AtmosphereShader.Use();

        for (int i = 0; i < 6; i++)
        {
            const glm::mat4& view_matrix = view_matrices[i];

            map.BindFace(i);
            m_AtmosphereShader.SetMatrix4("u_Projection", projection_matrix);
            m_AtmosphereShader.SetMatrix4("u_View", glm::mat4(glm::mat3(view_matrix)));
            m_AtmosphereShader.SetMatrix4("u_InvProjection", glm::inverse(projection_matrix));
            m_AtmosphereShader.SetMatrix4("u_InvView", glm::inverse(glm::mat4(glm::mat3(view_matrix))));
            m_AtmosphereShader.SetInteger("u_Skybox", 0);
            m_AtmosphereShader.SetFloat("u_Time", glfwGetTime());
            m_AtmosphereShader.SetVector3f("u_SunDirection", sun_direction);
            m_AtmosphereShader.SetInteger("u_NumSamples", steps);
            m_AtmosphereShader.SetInteger("u_NumLightSamples", lsteps);

            m_VAO.Bind();
            glDrawArrays(GL_TRIANGLES, 0, 6);
            m_VAO.Unbind();
        }

        glDepthMask(GL_TRUE);
    }

    void AtmosphereRenderer::Recompile()
    {
        m_AtmosphereShader.Destroy();
        m_AtmosphereShader.CreateShaderProgramFromFile("Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/AtmosphereFrag.glsl");
        m_AtmosphereShader.CompileShaders();
    }
}
