/* VoxelRT - A Voxel Raytracing Engine
Written by : Samuel Rasquinha

Contributors : 
UglySwedishFish
Telo
Moonsheep
Snurf

Resources : 
https://github.com/BrutPitt/glslSmartDeNoise/
https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
ScratchAPixel
*/

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
#include "Core/Player.h"
#include "Core/GLClasses/FramebufferRed.h"
#include "Core/ColorPassFBO.h"
#include "Core/GLClasses/Texture.h"
#include "Core/Renderer2D.h"

using namespace VoxelRT;
Player MainPlayer;
bool VSync = true;

float InitialTraceResolution = 0.75f;
float DiffuseTraceResolution = 0.25f;
float ShadowTraceResolution = 0.40;
float ReflectionTraceResolution = 0.25;
float SunTick = 50.0f;

bool FullyDynamicShadows = true;

glm::vec3 SunDirection;
glm::vec3 MoonDirection;

World* world = nullptr;
bool ModifiedWorld = false;
FPSCamera& MainCamera = MainPlayer.Camera;
VoxelRT::OrthographicCamera OCamera(0.0f, 800.0f, 0.0f, 600.0f);

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
		ImGui::SliderFloat("Initial Trace Resolution", &InitialTraceResolution, 0.1f, 1.25f);
		ImGui::SliderFloat("Diffuse Trace Resolution ", &DiffuseTraceResolution, 0.1f, 1.25f);
		ImGui::SliderFloat("Shadow Trace Resolution ", &ShadowTraceResolution, 0.1f, 1.25f);
		ImGui::SliderFloat("Reflection Trace Resolution ", &ReflectionTraceResolution, 0.1f, 0.8f);
		ImGui::SliderFloat("Sun Time ", &SunTick, 0.1f, 256.0f);
		ImGui::Checkbox("Fully Dynamic Shadows?", &FullyDynamicShadows);
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

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_F)
		{
			MainPlayer.Freefly = !MainPlayer.Freefly;
		}

		if (e.type == EventTypes::KeyPress && e.key == GLFW_KEY_Q)
		{
			world->ChangeCurrentlyHeldBlock();
		}

		if (e.type == EventTypes::WindowResize)
		{
			MainCamera.SetAspect((float)e.wx / (float)e.wy);
			OCamera.SetProjection(0.0f, e.wx, 0.0f, e.wy);
		}
	}

};

int main()
{
	RayTracerApp app;
	app.Initialize();
	BlockDatabase::Initialize();

	bool gen_type = 0;

	std::cout << "\nWhat type of world would you like to generate? (FLAT = 0, SIMPLEX FRACTAL = 1) : "; 
	std::cin >> gen_type;
	std::cout << "\n\n";

	// --
	world = new World();
	GenerateWorld(world, gen_type);
	world->Buffer();
	// --

	VoxelRT::Renderer2D RendererUI;

	GLClasses::VertexBuffer VBO;
	GLClasses::VertexArray VAO;
	GLClasses::Shader InitialTraceShader;
	GLClasses::Shader FinalShader;
	GLClasses::Shader DiffuseTraceShader;
	GLClasses::Shader DiffuseTemporalFilter;
	GLClasses::Shader DenoiseFilter;
	GLClasses::Shader ColorShader;
	GLClasses::Shader TemporalShader;
	GLClasses::Shader PostProcessingShader;
	GLClasses::Shader TemporalAAShader;
	GLClasses::Shader ShadowTraceShader;
	GLClasses::Shader ReflectionTraceShader;

	VoxelRT::InitialRTFBO InitialTraceFBO;
	GLClasses::Framebuffer DiffuseTraceFBO;
	GLClasses::Framebuffer DiffuseTemporalFBO1;
	GLClasses::Framebuffer DiffuseTemporalFBO2;
	GLClasses::Framebuffer DenoisedFBO;
	ColorPassFBO ColoredFBO;
	GLClasses::Framebuffer PostProcessingFBO;
	GLClasses::Framebuffer TAAFBO1;
	GLClasses::Framebuffer TAAFBO2;
	GLClasses::TextureArray BlueNoise;
	GLClasses::Framebuffer ShadowFBO;
	GLClasses::Framebuffer ReflectionTraceFBO;


	GLClasses::CubeTextureMap Skymap;
	GLClasses::CubeTextureMap SkymapLOWRES;
	glm::mat4 CurrentProjection, CurrentView;
	glm::mat4 PreviousProjection, PreviousView;
	glm::mat4 ShadowProjection, ShadowView;
	glm::mat4 ReflectionProjection, ReflectionView;
	glm::vec3 CurrentPosition, PreviousPosition;

	GLClasses::Texture Crosshair;
	Crosshair.CreateTexture("Res/Misc/crosshair.png", false);

	InitialTraceShader.CreateShaderProgramFromFile("Core/Shaders/InitialRayTraceVert.glsl", "Core/Shaders/InitialRayTraceFrag.glsl");
	InitialTraceShader.CompileShaders();
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
	TemporalAAShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalAA.glsl");
	TemporalAAShader.CompileShaders();
	ShadowTraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/ShadowRayTraceFrag.glsl");
	ShadowTraceShader.CompileShaders();
	ReflectionTraceShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionTraceFrag.glsl");
	ReflectionTraceShader.CompileShaders();

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

	BlueNoise.CreateArray({
		"Res/Misc/BL_0.png",
		"Res/Misc/BL_1.png",
		"Res/Misc/BL_2.png",
		"Res/Misc/BL_3.png"
		}, { 256, 256 }, true, false);

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

	MainCamera.SetPosition(glm::vec3(WORLD_SIZE_X / 2, 70, WORLD_SIZE_Z / 2));

	///// Set the block texture data uniforms    /////

	InitialTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec3 data;

		data.x = BlockDatabase::GetBlockTexture(i, BlockDatabase::BlockFaceType::Top);
		data.y = BlockDatabase::GetBlockNormalTexture(i, BlockDatabase::BlockFaceType::Top);
		data.z = BlockDatabase::GetBlockPBRTexture(i, BlockDatabase::BlockFaceType::Top);

		InitialTraceShader.SetVector3f(name.c_str(), data);
	}

	glUseProgram(0);

	DiffuseTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec3 data;

		data.x = BlockDatabase::GetBlockTexture(i, BlockDatabase::BlockFaceType::Top);
		data.y = BlockDatabase::GetBlockNormalTexture(i, BlockDatabase::BlockFaceType::Top);
		data.z = BlockDatabase::GetBlockPBRTexture(i, BlockDatabase::BlockFaceType::Top);

		DiffuseTraceShader.SetVector3f(name.c_str(), data);
	}

	glUseProgram(0);

	ReflectionTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec3 data;

		data.x = BlockDatabase::GetBlockTexture(i, BlockDatabase::BlockFaceType::Top);
		data.y = BlockDatabase::GetBlockNormalTexture(i, BlockDatabase::BlockFaceType::Top);
		data.z = BlockDatabase::GetBlockPBRTexture(i, BlockDatabase::BlockFaceType::Top);

		ReflectionTraceShader.SetVector3f(name.c_str(), data);
	}

	glUseProgram(0);

	/////                                     /////

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		// Tick the sun and moon
		float time_angle = SunTick * 2.0f;
		glm::mat4 sun_rotation_matrix;

		sun_rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(time_angle), glm::vec3(0.0f, 0.0f, 1.0f));
		SunDirection = glm::vec3(sun_rotation_matrix * glm::vec4(1.0f));
		MoonDirection = glm::vec3(-SunDirection.x, -SunDirection.y, SunDirection.z);


		glfwSwapInterval((int)VSync);

		// Resize the framebuffers
		DiffuseTraceFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		DiffuseTemporalFBO1.SetSize(app.GetWidth(), app.GetHeight());
		DiffuseTemporalFBO2.SetSize(app.GetWidth(), app.GetHeight());
		DenoisedFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		
		PostProcessingFBO.SetSize(app.GetWidth(), app.GetHeight());
		ColoredFBO.SetDimensions(app.GetWidth(), app.GetHeight());
		
		InitialTraceFBO.SetDimensions(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution));
		
		TAAFBO1.SetSize(app.GetWidth(), app.GetHeight());
		TAAFBO2.SetSize(app.GetWidth(), app.GetHeight());
		
		ShadowFBO.SetSize(app.GetWidth() * ShadowTraceResolution, app.GetHeight() * ShadowTraceResolution);
		
		ReflectionTraceFBO.SetSize(app.GetWidth() * ReflectionTraceResolution, app.GetHeight() * ReflectionTraceResolution);



		GLClasses::Framebuffer& TAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO1 : TAAFBO2;
		GLClasses::Framebuffer& PrevTAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO2 : TAAFBO1;

		GLClasses::Framebuffer& DiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO1 : DiffuseTemporalFBO2;
		GLClasses::Framebuffer& PrevDiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO2 : DiffuseTemporalFBO1;

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_F2) == GLFW_PRESS)
		{
			InitialTraceShader.Recompile();
			FinalShader.Recompile();
			DiffuseTraceShader.Recompile();
			DiffuseTemporalFilter.Recompile();
			DenoiseFilter.Recompile();
			ColorShader.Recompile();
			TemporalShader.Recompile();
			PostProcessingShader.Recompile();
			TemporalAAShader.Recompile();
			ShadowTraceShader.Recompile();
			ReflectionTraceShader.Recompile();

			///// Set the block texture data uniforms    /////

			InitialTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec3 data;

				data.x = BlockDatabase::GetBlockTexture(i, BlockDatabase::BlockFaceType::Top);
				data.y = BlockDatabase::GetBlockNormalTexture(i, BlockDatabase::BlockFaceType::Top);
				data.z = BlockDatabase::GetBlockPBRTexture(i, BlockDatabase::BlockFaceType::Top);

				InitialTraceShader.SetVector3f(name.c_str(), data);
			}

			glUseProgram(0);

			DiffuseTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec3 data;

				data.x = BlockDatabase::GetBlockTexture(i, BlockDatabase::BlockFaceType::Top);
				data.y = BlockDatabase::GetBlockNormalTexture(i, BlockDatabase::BlockFaceType::Top);
				data.z = BlockDatabase::GetBlockPBRTexture(i, BlockDatabase::BlockFaceType::Top);

				DiffuseTraceShader.SetVector3f(name.c_str(), data);
			}

			glUseProgram(0);

			ReflectionTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec3 data;

				data.x = BlockDatabase::GetBlockTexture(i, BlockDatabase::BlockFaceType::Top);
				data.y = BlockDatabase::GetBlockNormalTexture(i, BlockDatabase::BlockFaceType::Top);
				data.z = BlockDatabase::GetBlockPBRTexture(i, BlockDatabase::BlockFaceType::Top);

				ReflectionTraceShader.SetVector3f(name.c_str(), data);
			}

			glUseProgram(0);

			/////                                     /////

			Logger::Log("Recompiled!");
		}

		MainPlayer.OnUpdate(app.GetWindow(), world);
		app.OnUpdate();

		// ---------------------

		glm::mat4 TempView = PreviousView;

		PreviousProjection = CurrentProjection;
		PreviousView = CurrentView;
		PreviousPosition = CurrentPosition;
		CurrentProjection = MainCamera.GetProjectionMatrix();
		CurrentView = MainCamera.GetViewMatrix();
		CurrentPosition = MainCamera.GetPosition();

		glm::mat4 inv_view = glm::inverse(MainCamera.GetViewMatrix());
		glm::mat4 inv_projection = glm::inverse(MainCamera.GetProjectionMatrix());
		bool PlayerMoved = TempView != MainCamera.GetViewMatrix();

		if (TempView != MainCamera.GetViewMatrix() || ModifiedWorld || app.GetCurrentFrame() % 4 == 0)
		{
			InitialTraceFBO.Bind();
			glDisable(GL_CULL_FACE);
			glDisable(GL_DEPTH_TEST);

			InitialTraceShader.Use();

			InitialTraceShader.SetMatrix4("u_InverseView", inv_view);
			InitialTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
			InitialTraceShader.SetInteger("u_VoxelDataTexture", 0);
			InitialTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			InitialTraceShader.SetInteger("u_VertCurrentFrame", app.GetCurrentFrame());
			InitialTraceShader.SetVector2f("u_Dimensions", glm::vec2(InitialTraceFBO.GetWidth(), InitialTraceFBO.GetHeight()));
			InitialTraceShader.SetVector2f("u_VertDimensions", glm::vec2(app.GetWidth(), app.GetHeight()));

			/// TEMPORARY ///
			InitialTraceShader.SetInteger("u_GrassBlockProps[0]", BlockDatabase::GetBlockID("Grass"));
			InitialTraceShader.SetInteger("u_GrassBlockProps[1]", BlockDatabase::GetBlockTexture("Grass", BlockDatabase::BlockFaceType::Top));
			InitialTraceShader.SetInteger("u_GrassBlockProps[2]", BlockDatabase::GetBlockNormalTexture("Grass", BlockDatabase::BlockFaceType::Top));
			InitialTraceShader.SetInteger("u_GrassBlockProps[3]", BlockDatabase::GetBlockPBRTexture("Grass", BlockDatabase::BlockFaceType::Top));
			InitialTraceShader.SetInteger("u_GrassBlockProps[4]", BlockDatabase::GetBlockTexture("Grass", BlockDatabase::BlockFaceType::Front));
			InitialTraceShader.SetInteger("u_GrassBlockProps[5]", BlockDatabase::GetBlockNormalTexture("Grass", BlockDatabase::BlockFaceType::Front));
			InitialTraceShader.SetInteger("u_GrassBlockProps[6]", BlockDatabase::GetBlockPBRTexture("Grass", BlockDatabase::BlockFaceType::Front));
			InitialTraceShader.SetInteger("u_GrassBlockProps[7]", BlockDatabase::GetBlockTexture("Grass", BlockDatabase::BlockFaceType::Bottom));
			InitialTraceShader.SetInteger("u_GrassBlockProps[8]", BlockDatabase::GetBlockNormalTexture("Grass", BlockDatabase::BlockFaceType::Bottom));
			InitialTraceShader.SetInteger("u_GrassBlockProps[9]", BlockDatabase::GetBlockPBRTexture("Grass", BlockDatabase::BlockFaceType::Bottom));

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
		glClear(GL_COLOR_BUFFER_BIT);

		DiffuseTraceShader.Use();

		DiffuseTraceShader.SetInteger("u_VoxelData", 0);
		DiffuseTraceShader.SetInteger("u_PositionTexture", 1);
		DiffuseTraceShader.SetInteger("u_NormalTexture", 2);
		DiffuseTraceShader.SetInteger("u_Skymap", 3);
		DiffuseTraceShader.SetInteger("u_DataTexture", 5);
		DiffuseTraceShader.SetInteger("u_BlockNormalTextures", 4);
		DiffuseTraceShader.SetInteger("u_BlockAlbedoTextures", 6);
		DiffuseTraceShader.SetInteger("u_BlueNoiseTextures", 7);
		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
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

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetDataTexture());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlueNoise.GetTextureArray());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTraceFBO.Unbind();


		// ---- Denoise the diffuse ----

		DenoisedFBO.Bind();
		DenoiseFilter.Use();

		DenoiseFilter.SetInteger("u_Texture", 0);
		DenoiseFilter.SetInteger("u_NormalTexture", 1);
		DenoiseFilter.SetInteger("u_InitialTracePositionTexture", 2);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetNormalTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTraceFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DenoisedFBO.Unbind();

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
		DiffuseTemporalFilter.SetBool("u_CameraMoved", PreviousPosition != MainCamera.GetPosition());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DenoisedFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTemporalFBO.Unbind();

		// ---- SHADOW TRACE ----

		bool UpdateShadows = FullyDynamicShadows ? true : (app.GetCurrentFrame() % 4 == 0);

		if (ModifiedWorld || UpdateShadows)
		{
			ShadowProjection = CurrentProjection;
			ShadowView = CurrentView;

			ShadowFBO.Bind();
			ShadowTraceShader.Use();

			ShadowTraceShader.SetInteger("u_PositionTexture", 0);
			ShadowTraceShader.SetInteger("u_VoxelData", 1);
			ShadowTraceShader.SetVector3f("u_SunDirection", SunDirection);
			ShadowTraceShader.SetVector3f("u_MoonDirection", MoonDirection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ShadowFBO.Unbind();
		}

		// ---- COLOR PASS ----

		ReflectionProjection = PreviousProjection;
		ReflectionView = PreviousView;

		ColoredFBO.Bind();

		ColorShader.Use();
		ColorShader.SetInteger("u_DiffuseTexture", 0);
		ColorShader.SetInteger("u_NormalTexture", 1);
		ColorShader.SetInteger("u_InitialTracePositionTexture", 2);
		ColorShader.SetInteger("u_DataTexture", 3);
		ColorShader.SetInteger("u_BlockAlbedoTextures", 4);
		ColorShader.SetInteger("u_BlockNormalTextures", 5);
		ColorShader.SetInteger("u_BlockPBRTextures", 6);
		ColorShader.SetInteger("u_Skybox", 7);
		ColorShader.SetInteger("u_ShadowTexture", 8);
		ColorShader.SetInteger("u_BlueNoiseTextures", 9);
		ColorShader.SetInteger("u_ReflectionTraceTexture", 10);
		ColorShader.SetMatrix4("u_InverseView", inv_view);
		ColorShader.SetMatrix4("u_InverseProjection", inv_projection);
		ColorShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
		ColorShader.SetMatrix4("u_ShadowView", ShadowView);
		ColorShader.SetMatrix4("u_ReflectionProjection", ReflectionProjection);
		ColorShader.SetMatrix4("u_ReflectionView", ReflectionView);
		ColorShader.SetVector2f("u_InitialTraceResolution", glm::vec2(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution)));
		ColorShader.SetVector3f("u_SunDirection", SunDirection);
		ColorShader.SetVector3f("u_MoonDirection", MoonDirection);
		ColorShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetNormalTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetDataTexture());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetPBRTextureArray());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, ShadowFBO.GetTexture());

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlueNoise.GetTextureArray());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		ColoredFBO.Unbind();

		// ----- REFLECTION TRACE -----

		{
			ReflectionTraceFBO.Bind();
			ReflectionTraceShader.Use();

			ReflectionTraceShader.SetInteger("u_PositionTexture", 0);
			ReflectionTraceShader.SetInteger("u_NormalTexture", 1);
			ReflectionTraceShader.SetInteger("u_PBRTexture", 2);
			ReflectionTraceShader.SetInteger("u_VoxelData", 3);
			ReflectionTraceShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
			ReflectionTraceShader.SetInteger("u_BlockNormalTextures", 4);
			ReflectionTraceShader.SetInteger("u_BlockAlbedoTextures", 5);
			ReflectionTraceShader.SetInteger("u_Skymap", 6);
			ReflectionTraceShader.SetInteger("u_InitialTraceNormalTexture", 7);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetNormalTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetPBRTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE6);
			glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapLOWRES.GetID());

			glActiveTexture(GL_TEXTURE7);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetNormalTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionTraceFBO.Unbind();
		}

		// ----- TAA ----- //

		TAAFBO.Bind();

		TemporalAAShader.Use();
		TemporalAAShader.SetInteger("u_CurrentColorTexture", 0);
		TemporalAAShader.SetInteger("u_PositionTexture", 1);
		TemporalAAShader.SetInteger("u_PreviousColorTexture", 2);

		TemporalAAShader.SetMatrix4("u_PrevProjection", PreviousProjection);
		TemporalAAShader.SetMatrix4("u_PrevView", PreviousView);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetColorTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO.GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevTAAFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		// ---- POST PROCESSING ----

		PostProcessingFBO.Bind();
		PostProcessingShader.Use();

		PostProcessingShader.SetMatrix4("u_InverseView", inv_view);
		PostProcessingShader.SetMatrix4("u_InverseProjection", inv_projection);
		PostProcessingShader.SetInteger("u_FramebufferTexture", 0);
		PostProcessingShader.SetInteger("u_PositionTexture", 1);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TAAFBO.GetTexture());

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

		RendererUI.RenderQuad(glm::vec2(floor((float)app.GetWidth() / 2.0f), floor((float)app.GetHeight() / 2.0f)), &Crosshair, &OCamera);

		// Finish Frame
		app.FinishFrame();
		GLClasses::DisplayFrameRate(app.GetWindow(), "Voxel RT");

		ModifiedWorld = false;
	}

	delete world;
}