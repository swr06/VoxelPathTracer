#include "AtmosphereRenderer.h"

#include <glm/glm.hpp>

namespace VoxelRT
{
    AtmosphereRenderer::AtmosphereRenderer() : m_VBO(GL_ARRAY_BUFFER)
	{
        m_AtmosphereShader.CreateShaderProgramFromFile("Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/AtmosphereFrag.glsl");
        m_AtmosphereShader.CompileShaders();
        m_AtmosphereNebulaCombine.CreateShaderProgramFromFile("Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/CubemapCombine.glsl");
        m_AtmosphereNebulaCombine.CompileShaders();

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

    void AtmosphereRenderer::DownsampleAtmosphere(AtmosphereRenderMap& map, AtmosphereRenderMap& atmos_map, const glm::mat4& Matrix, GLuint Texture, float SunVisibility, float strength)
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

        m_AtmosphereNebulaCombine.Use();

        for (int i = 0; i < 6; i++)
        {
            const glm::mat4& view_matrix = view_matrices[i];

            map.Bind();
            map.BindFace(i);


            m_AtmosphereNebulaCombine.SetMatrix4("u_Projection", projection_matrix);
            m_AtmosphereNebulaCombine.SetMatrix4("u_View", glm::mat4(glm::mat3(view_matrix)));
            m_AtmosphereNebulaCombine.SetMatrix4("u_InvProjection", glm::inverse(projection_matrix));
            m_AtmosphereNebulaCombine.SetMatrix4("u_InvView", glm::inverse(glm::mat4(glm::mat3(view_matrix))));

            m_AtmosphereNebulaCombine.SetInteger("u_Skybox", 0);
            m_AtmosphereNebulaCombine.SetInteger("u_NebulaLowRes", 1);
            m_AtmosphereNebulaCombine.SetMatrix4("u_RotationMatrix", Matrix);
            m_AtmosphereNebulaCombine.SetFloat("u_SunVisibility", SunVisibility);
            m_AtmosphereNebulaCombine.SetFloat("u_SunVisib", SunVisibility);
            m_AtmosphereNebulaCombine.SetFloat("u_NebulaStrength", strength);
            
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_CUBE_MAP, atmos_map.GetTexture());

            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_CUBE_MAP, Texture);

            m_VAO.Bind();
            glDrawArrays(GL_TRIANGLES, 0, 6);
            m_VAO.Unbind();
        }

    }


    void AtmosphereRenderer::Recompile()
    {
        m_AtmosphereShader.Destroy();
        m_AtmosphereShader.CreateShaderProgramFromFile("Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/AtmosphereFrag.glsl");
        m_AtmosphereShader.CompileShaders();

        m_AtmosphereNebulaCombine.Destroy();
        m_AtmosphereNebulaCombine.CreateShaderProgramFromFile("Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/CubemapCombine.glsl");
        m_AtmosphereNebulaCombine.CompileShaders();
    }
}
