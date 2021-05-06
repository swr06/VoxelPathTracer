#include <iostream>

#include "Core/Application/Application.h"
#include "Core/GLClasses/VertexArray.h"
#include "Core/GLClasses/VertexBuffer.h"
#include "Core/GLClasses/IndexBuffer.h"
#include "Core/GLClasses/Shader.h"
#include "Core/FpsCamera.h"
#include "Core/GLClasses/Fps.h"
#include "Core/World.h"
#include "Core/GLClasses/Framebuffer.h"
#include "Core/WorldGenerator.h"
#include "Core/InitialTraceFBO.h"
#include "Core/GLClasses/CubeTextureMap.h"
#include "Core/BlockDatabase.h"

using namespace VoxelRT;
FPSCamera MainCamera(90.0f, (float)800.0f / (float)600.0f);
bool VSync = true;

float InitialTraceResolution = 1.0f;
float DiffuseTraceResolution = 0.5f;

class RayTracerApp : public Application
{
public:

	RayTracerApp()
	{
		m_Width = 800;
		m_Height = 600;
	}

	void OnUserCreate(double ts) override
	{

	}

	void OnUserUpdate(double ts) override
	{

	}

	void OnImguiRender(double ts) override
	{
		ImGui::Text("Player Position : %f, %f, %f", MainCamera.GetPosition().x, MainCamera.GetPosition().y, MainCamera.GetPosition().z);
		ImGui::Text("Camera Front : %f, %f, %f", MainCamera.GetFront().x, MainCamera.GetFront().y, MainCamera.GetFront().z);
	}

	void OnEvent(Event e) override
	{
		if (e.type == EventTypes::MouseMove && GetCursorLocked())
		{
			MainCamera.UpdateOnMouseMovement(GetCursorX(), GetCursorY());
		}

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_F1)
		{
			this->SetCursorLocked(!this->GetCursorLocked());
		}

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_V)
		{
			VSync = !VSync;
		}

		if (e.type == EventTypes::WindowResize)
		{
			MainCamera.SetAspect((float)e.wx / (float)e.wy);
		}
	}

};

int main()
{
	RayTracerApp app;
	app.Initialize();
	BlockDatabase::Initialize();

	GLClasses::VertexBuffer VBO;
	GLClasses::VertexArray VAO;
	GLClasses::Shader RaytraceShader;
	GLClasses::Shader FinalShader;
	GLClasses::Shader DiffuseTraceShader;
	GLClasses::Shader TemporalFilter;
	GLClasses::Shader DenoiseFilter;
	GLClasses::Shader ColorShader;

	VoxelRT::InitialRTFBO InitialTraceFBO;
	GLClasses::Framebuffer DiffuseTraceFBO;
	GLClasses::Framebuffer TemporalFBO1;
	GLClasses::Framebuffer TemporalFBO2;
	GLClasses::Framebuffer DenoisedFBO;
	GLClasses::Framebuffer ColoredFBO;

	GLClasses::CubeTextureMap Skymap;
	glm::mat4 CurrentProjection, CurrentView;
	glm::mat4 PreviousProjection, PreviousView;

	RaytraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/InitialRayTraceFrag.glsl");
	RaytraceShader.CompileShaders();
	FinalShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/FBOFrag.glsl");
	FinalShader.CompileShaders();
	DiffuseTraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/DiffuseRayTraceFrag.glsl");
	DiffuseTraceShader.CompileShaders();
	TemporalFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalFilter.glsl");
	TemporalFilter.CompileShaders();
	DenoiseFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SmartDenoise.glsl");
	DenoiseFilter.CompileShaders();
	ColorShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/ColorPassFrag.glsl");
	ColorShader.CompileShaders();

	Skymap.CreateCubeTextureMap(
		{
		"Res/right.bmp",
		"Res/left.bmp",
		"Res/top.bmp",
		"Res/bottom.bmp",
		"Res/front.bmp",
		"Res/back.bmp"
		}, false
	);

	float Vertices[] =
	{
		-1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
		 1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
		 1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
	};

	VAO.Bind();
	VBO.Bind();
	VBO.BufferData(sizeof(Vertices), Vertices, GL_STATIC_DRAW);
	VBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
	VBO.VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
	VAO.Unbind();

	// ----

	World* world = new World();
	GenerateWorld(world);
	world->Buffer();

	// ----

	app.SetCursorLocked(true);

	glDisable(GL_BLEND);

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		GLClasses::Framebuffer& TemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? TemporalFBO1 : TemporalFBO2;
		GLClasses::Framebuffer& PrevTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? TemporalFBO2 : TemporalFBO1;

		glfwSwapInterval((int)VSync);
		InitialTraceFBO.SetDimensions(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution));
		DiffuseTraceFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		TemporalFBO.SetSize(app.GetWidth(), app.GetHeight());
		PrevTemporalFBO.SetSize(app.GetWidth(), app.GetHeight());
		ColoredFBO.SetSize(app.GetWidth(), app.GetHeight());
		DenoisedFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);

		float camera_speed = 0.125f;

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_W) == GLFW_PRESS)
		{
			// Take the cross product of the MainCamera's right and up.
			MainCamera.ChangePosition(MainCamera.GetFront() * camera_speed);
		}

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_S) == GLFW_PRESS)
		{
			MainCamera.ChangePosition(-MainCamera.GetFront() * camera_speed);
		}

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_A) == GLFW_PRESS)
		{
			MainCamera.ChangePosition(-(MainCamera.GetRight() * camera_speed));
		}

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_D) == GLFW_PRESS)
		{
			MainCamera.ChangePosition(MainCamera.GetRight() * camera_speed);
		}

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_SPACE) == GLFW_PRESS)
		{
			MainCamera.ChangePosition(MainCamera.GetUp() * camera_speed);
		}

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		{
			MainCamera.ChangePosition(-(MainCamera.GetUp() * camera_speed));
		}

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_F2) == GLFW_PRESS)
		{
			RaytraceShader.Recompile();
			DiffuseTraceShader.Recompile();
			FinalShader.Recompile();
			TemporalFilter.Recompile();
			DenoiseFilter.Recompile();
			ColorShader.Recompile();
			Logger::Log("Recompiled!");
		}

		MainCamera.OnUpdate();
		MainCamera.Refresh();

		app.OnUpdate();

		// ---------------------

		PreviousProjection = CurrentProjection;
		PreviousView = CurrentView;
		CurrentProjection = MainCamera.GetProjectionMatrix();
		CurrentView = MainCamera.GetViewMatrix();

		InitialTraceFBO.Bind();
		glClear(GL_COLOR_BUFFER_BIT);
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		glm::mat4 inv_view = glm::inverse(MainCamera.GetViewMatrix());
		glm::mat4 inv_projection = glm::inverse(MainCamera.GetProjectionMatrix());

		RaytraceShader.Use();

		RaytraceShader.SetMatrix4("u_InverseView", inv_view);
		RaytraceShader.SetMatrix4("u_InverseProjection", inv_projection);
		RaytraceShader.SetInteger("u_VoxelDataTexture", 0);
		RaytraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		// Diffuse tracing

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		DiffuseTraceFBO.Bind();
		DiffuseTraceShader.Use();

		DiffuseTraceShader.SetInteger("u_VoxelData", 0);
		DiffuseTraceShader.SetInteger("u_PositionTexture", 1);
		DiffuseTraceShader.SetInteger("u_NormalTexture", 2);
		DiffuseTraceShader.SetInteger("u_Skymap", 3);
		DiffuseTraceShader.SetMatrix4("u_InverseView", inv_view);
		DiffuseTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
		DiffuseTraceShader.SetVector2f("u_Dimensions", glm::vec2(DiffuseTraceFBO.GetWidth(), DiffuseTraceFBO.GetHeight()));
		DiffuseTraceShader.SetFloat("u_Time", glfwGetTime() * 1.2f);
		RaytraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetNormalTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTraceFBO.Unbind();

		///// Temporal filtering /////

		TemporalFBO.Bind();
		TemporalFilter.Use();

		TemporalFilter.SetInteger("u_CurrentColorTexture", 0);
		TemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
		TemporalFilter.SetInteger("u_PreviousColorTexture", 2);

		TemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
		TemporalFilter.SetMatrix4("u_View", CurrentView);
		TemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
		TemporalFilter.SetMatrix4("u_PrevView", PreviousView);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTraceFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevTemporalFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		TemporalFBO.Unbind();

		// ---- Smart Denoise ----

		DenoisedFBO.Bind();
		DenoiseFilter.Use();
		
		DenoiseFilter.SetInteger("u_Texture", 0);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TemporalFBO.GetTexture());
		
		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();
		
		DenoisedFBO.Unbind();

		// ---- COLOR PASS ----

		ColoredFBO.Bind();

		ColorShader.Use();
		ColorShader.SetInteger("u_DiffuseTexture", 0);
		ColorShader.SetInteger("u_NormalTexture", 1);
		ColorShader.SetInteger("u_InitialTracePositionTexture", 2);
		ColorShader.SetInteger("u_DataTexture", 3);
		ColorShader.SetInteger("u_BlockAlbedoTextures", 4);
		ColorShader.SetInteger("u_Skybox", 5);
		ColorShader.SetMatrix4("u_InverseView", inv_view);
		ColorShader.SetMatrix4("u_InverseProjection", inv_projection);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DenoisedFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetNormalTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetDataTexture());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		ColoredFBO.Unbind();


		// ---- FINAL ----

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, app.GetWidth(), app.GetHeight());

		FinalShader.Use();
		FinalShader.SetInteger("u_FramebufferTexture", 0);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		// Finish Frame
		app.FinishFrame();
		GLClasses::DisplayFrameRate(app.GetWindow(), "Voxel RT");
	}
}