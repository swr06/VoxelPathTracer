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

float InitialTraceResolution = 0.5f;
float DiffuseTraceResolution = 0.5f;
World* world = nullptr;
bool ModifiedWorld = false;

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

		if (e.type == EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_LEFT && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(false, MainCamera.GetPosition(), MainCamera.GetFront());
				ModifiedWorld = true;
			}
		}

		if (e.type == EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_RIGHT && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(true, MainCamera.GetPosition(), MainCamera.GetFront());
				ModifiedWorld = true;
			}
		}

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_F1)
		{
			this->SetCursorLocked(!this->GetCursorLocked());
		}

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_V)
		{
			VSync = !VSync;
		}

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_Q)
		{
			world->ChangeCurrentlyHeldBlock();
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

	// --
	world = new World();
	GenerateWorld(world);
	world->Buffer();

	// --


	GLClasses::VertexBuffer VBO;
	GLClasses::VertexArray VAO;
	GLClasses::Shader RaytraceShader;
	GLClasses::Shader FinalShader;
	GLClasses::Shader DiffuseTraceShader;
	GLClasses::Shader DiffuseTemporalFilter;
	GLClasses::Shader DenoiseFilter;
	GLClasses::Shader ColorShader;
	GLClasses::Shader TemporalShader;
	GLClasses::Shader PostProcessingShader;

	VoxelRT::InitialRTFBO InitialTraceFBO;
	GLClasses::Framebuffer DiffuseTraceFBO;
	GLClasses::Framebuffer DiffuseTemporalFBO1;
	GLClasses::Framebuffer DiffuseTemporalFBO2;
	GLClasses::Framebuffer DenoisedFBO;
	GLClasses::Framebuffer ColoredFBO;
	GLClasses::Framebuffer PostProcessingFBO;

	GLClasses::CubeTextureMap Skymap;
	GLClasses::CubeTextureMap SkymapLOWRES;
	glm::mat4 CurrentProjection, CurrentView;
	glm::mat4 PreviousProjection, PreviousView;

	RaytraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/InitialRayTraceFrag.glsl");
	RaytraceShader.CompileShaders();
	FinalShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/FBOFrag.glsl");
	FinalShader.CompileShaders();
	DiffuseTraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/DiffuseRayTraceFrag.glsl");
	DiffuseTraceShader.CompileShaders();
	DiffuseTemporalFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/DiffuseTemporalFilter.glsl");
	DiffuseTemporalFilter.CompileShaders();
	DenoiseFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SmartDenoise.glsl");
	DenoiseFilter.CompileShaders();
	ColorShader.CreateShaderProgramFromFile("Core/Shaders/ColorPassVert.glsl", "Core/Shaders/ColorPassFrag.glsl");
	ColorShader.CompileShaders();
	TemporalShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/ColorTemporalFilter.glsl");
	TemporalShader.CompileShaders();
	PostProcessingShader.CreateShaderProgramFromFile("Core/Shaders/PostProcessingVert.glsl", "Core/Shaders/PostProcessingFrag.glsl");
	PostProcessingShader.CompileShaders();

	Skymap.CreateCubeTextureMap(
		{
		"Res/right.bmp",
		"Res/left.bmp",
		"Res/top.bmp",
		"Res/bottom.bmp",
		"Res/front.bmp",
		"Res/back.bmp"
		}, true
	);

	SkymapLOWRES.CreateCubeTextureMap(
		{
		"Res/skymap_lowres/right.bmp",
		"Res/skymap_lowres/left.bmp",
		"Res/skymap_lowres/top.bmp",
		"Res/skymap_lowres/bottom.bmp",
		"Res/skymap_lowres/front.bmp",
		"Res/skymap_lowres/back.bmp"
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

	app.SetCursorLocked(true);

	glDisable(GL_BLEND);

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		glfwSwapInterval((int)VSync);
		
		// Resize the framebuffers
		DiffuseTraceFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		DiffuseTemporalFBO1.SetSize(app.GetWidth(), app.GetHeight());
		DiffuseTemporalFBO2.SetSize(app.GetWidth(), app.GetHeight());
		PostProcessingFBO.SetSize(app.GetWidth(), app.GetHeight());
		DenoisedFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		InitialTraceFBO.SetDimensions(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution));
		ColoredFBO.SetSize(app.GetWidth(), app.GetHeight());

		GLClasses::Framebuffer& DiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO1 : DiffuseTemporalFBO2;
		GLClasses::Framebuffer& PrevDiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO2 : DiffuseTemporalFBO1;

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
			FinalShader.Recompile();
			DiffuseTraceShader.Recompile();
			DiffuseTemporalFilter.Recompile();
			DenoiseFilter.Recompile();
			ColorShader.Recompile();
			TemporalShader.Recompile();
			Logger::Log("Recompiled!");
		}

		MainCamera.Refresh();

		app.OnUpdate();

		// ---------------------

		glm::mat4 TempView = PreviousView;
		PreviousProjection = CurrentProjection;
		PreviousView = CurrentView;
		CurrentProjection = MainCamera.GetProjectionMatrix();
		CurrentView = MainCamera.GetViewMatrix();

		glm::mat4 inv_view = glm::inverse(MainCamera.GetViewMatrix());
		glm::mat4 inv_projection = glm::inverse(MainCamera.GetProjectionMatrix());

		if (TempView != MainCamera.GetViewMatrix() || ModifiedWorld)
		{
			InitialTraceFBO.Bind();
			glDisable(GL_CULL_FACE);
			glDisable(GL_DEPTH_TEST);

			RaytraceShader.Use();

			RaytraceShader.SetMatrix4("u_InverseView", inv_view);
			RaytraceShader.SetMatrix4("u_InverseProjection", inv_projection);
			RaytraceShader.SetInteger("u_VoxelDataTexture", 0);
			RaytraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			RaytraceShader.SetInteger("u_VertCurrentFrame", app.GetCurrentFrame());
			RaytraceShader.SetVector2f("u_Dimensions", glm::vec2(InitialTraceFBO.GetWidth(), InitialTraceFBO.GetHeight()));
			RaytraceShader.SetVector2f("u_VertDimensions", glm::vec2(app.GetWidth(), app.GetHeight()));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			InitialTraceFBO.Unbind();
		}

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
		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetNormalTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapLOWRES.GetID());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTraceFBO.Unbind();

		///// Temporal filter the calculated diffuse /////

		DiffuseTemporalFBO.Bind();
		DiffuseTemporalFilter.Use();

		DiffuseTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
		DiffuseTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
		DiffuseTemporalFilter.SetInteger("u_PreviousColorTexture", 2);

		DiffuseTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
		DiffuseTemporalFilter.SetMatrix4("u_View", CurrentView);
		DiffuseTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
		DiffuseTemporalFilter.SetMatrix4("u_PrevView", PreviousView);
		
		DiffuseTemporalFilter.SetBool("u_WorldModified", ModifiedWorld);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTraceFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTemporalFBO.Unbind();

		// ---- Denoise the diffuse ----

		DenoisedFBO.Bind();
		DenoiseFilter.Use();
		
		DenoiseFilter.SetInteger("u_Texture", 0);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture());
		
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
		ColorShader.SetVector2f("u_InitialTraceResolution", glm::vec2(floor(app.GetWidth()* InitialTraceResolution), floor(app.GetHeight()* InitialTraceResolution)));

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

		// ---- POST PROCESSING ----

		PostProcessingFBO.Bind();
		PostProcessingShader.Use();

		PostProcessingShader.SetMatrix4("u_InverseView", inv_view);
		PostProcessingShader.SetMatrix4("u_InverseProjection", inv_projection);
		PostProcessingShader.SetInteger("u_FramebufferTexture", 0);
		PostProcessingShader.SetInteger("u_PositionTexture", 1);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		PostProcessingFBO.Unbind();

		// ---- FINAL ----

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, app.GetWidth(), app.GetHeight());

		FinalShader.Use();
		FinalShader.SetInteger("u_FramebufferTexture", 0);
		FinalShader.SetInteger("u_InitialTracePositionTexture", 1);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, PostProcessingFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		// Finish Frame
		app.FinishFrame();
		GLClasses::DisplayFrameRate(app.GetWindow(), "Voxel RT");

		ModifiedWorld = false;
	}

	delete world;
}