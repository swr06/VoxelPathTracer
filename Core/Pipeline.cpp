#include "Pipeline.h"

static VoxelRT::Player MainPlayer;
static bool VSync = false;

static float InitialTraceResolution = 0.75f;
static float DiffuseTraceResolution = 0.125f;

static float ShadowTraceResolution = 0.50;
static float ReflectionTraceResolution = 0.3;
static float SSAOResolution = 0.35f;

static float VolumetricResolution = 0.5f;

static float SunTick = 50.0f;
static float DiffuseLightIntensity = 80.0f;
static float LensFlareIntensity = 0.075f;

static bool Bloom = true;

static bool GodRays = false;
static bool FakeGodRays = false;

static bool LensFlare = false;
static bool SSAO = false;

static bool FullyDynamicShadows = false;

static int GodRaysStepCount = 12;

static bool AutoExposure = false;

static glm::vec3 SunDirection;
static glm::vec3 MoonDirection;

static VoxelRT::World* world = nullptr;
static bool ModifiedWorld = false;
static VoxelRT::FPSCamera& MainCamera = MainPlayer.Camera;
static VoxelRT::OrthographicCamera OCamera(0.0f, 800.0f, 0.0f, 600.0f);

class RayTracerApp : public VoxelRT::Application
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
		ImGui::SliderFloat("SSAO Render Resolution ", &SSAOResolution, 0.1f, 0.9f);
		ImGui::SliderFloat("Diffuse Light Intensity ", &DiffuseLightIntensity, -10.0, 80.0f);
		ImGui::SliderFloat("Sun Time ", &SunTick, 0.1f, 256.0f);
		ImGui::SliderFloat("Lens Flare Intensity ", &LensFlareIntensity, 0.05f, 1.25f);
		ImGui::SliderInt("God ray raymarch step count", &GodRaysStepCount, 8, 64);
		ImGui::Checkbox("Lens Flare?", &LensFlare);
		ImGui::Checkbox("Fully Dynamic Shadows? (Fixes shadow artifacts)", &FullyDynamicShadows);
		ImGui::Checkbox("(Implementation - 1) God Rays? (Slower)", &GodRays);
		ImGui::Checkbox("(Implementation - 2) God Rays? (faster, more crisp, Adjust the step count in the menu)", &FakeGodRays);
		ImGui::Checkbox("Screen Space Ambient Occlusion?", &SSAO);
		ImGui::Checkbox("Bloom (Expensive!) ?", &Bloom);
		ImGui::Checkbox("Auto Exposure (Very very WIP!) ?", &AutoExposure);
	}

	void OnEvent(VoxelRT::Event e) override
	{
		if (e.type == VoxelRT::EventTypes::MouseMove && GetCursorLocked())
		{
			MainCamera.UpdateOnMouseMovement(GetCursorX(), GetCursorY());
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_LEFT && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(false, MainCamera.GetPosition(), MainCamera.GetFront());
				ModifiedWorld = true;
			}
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_RIGHT && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(true, MainCamera.GetPosition(), MainCamera.GetFront());
				ModifiedWorld = true;
			}
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_F1)
		{
			this->SetCursorLocked(!this->GetCursorLocked());
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_V)
		{
			VSync = !VSync;
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_F)
		{
			MainPlayer.Freefly = !MainPlayer.Freefly;
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_Q)
		{
			world->ChangeCurrentlyHeldBlock(true);
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_E)
		{
			world->ChangeCurrentlyHeldBlock(false);
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_ESCAPE)
		{
			VoxelRT::SaveWorld(world, world->m_Name);
			delete world;
			exit(0);
		}

		if (e.type == VoxelRT::EventTypes::WindowResize)
		{
			MainCamera.SetAspect((float)e.wx / (float)e.wy);
			OCamera.SetProjection(0.0f, e.wx, 0.0f, e.wy);
		}
	}

};

void VoxelRT::MainPipeline::StartPipeline()
{
	RayTracerApp app;
	app.Initialize();
	VoxelRT::BlockDatabase::Initialize();

	bool gen_type = 0;

	std::string world_name;

	world = new VoxelRT::World();

	do {
		std::cout << "\nEnter the name of your world : ";
		std::cin >> world_name;
	} while (!VoxelRT::FilenameValid(world_name));

	world->m_Name = world_name;

	if (!LoadWorld(world, world_name))
	{
		std::cout << "\nWhat type of world would you like to generate? (FLAT = 0, SIMPLEX FRACTAL = 1) : ";
		std::cin >> gen_type;
		std::cout << "\n\n";

		GenerateWorld(world, gen_type);
	}

	world->Buffer();



	VoxelRT::Renderer2D RendererUI;

	GLClasses::VertexBuffer VBO;
	GLClasses::VertexArray VAO;
	GLClasses::Shader InitialTraceShader;
	GLClasses::Shader FinalShader;
	GLClasses::Shader DiffuseTraceShader;
	GLClasses::Shader TemporalFilter;
	GLClasses::Shader DenoiseFilter;
	GLClasses::Shader ColorShader;
	GLClasses::Shader PostProcessingShader;
	GLClasses::Shader TemporalAAShader;
	GLClasses::Shader ShadowTraceShader;
	GLClasses::Shader ReflectionTraceShader;
	GLClasses::Shader SSAOShader;
	GLClasses::Shader SSAO_Blur;
	GLClasses::Shader SimpleDownsample;
	GLClasses::Shader LumaAverager;
	GLClasses::Shader VolumetricScattering;
	GLClasses::Shader BilateralBlur;
	GLClasses::Shader ReflectionDenoiser;

	VoxelRT::InitialRTFBO InitialTraceFBO_1;
	VoxelRT::InitialRTFBO InitialTraceFBO_2;
	GLClasses::Framebuffer DiffuseTraceFBO;
	GLClasses::Framebuffer DiffuseTemporalFBO1;
	GLClasses::Framebuffer DiffuseTemporalFBO2;
	GLClasses::Framebuffer DiffuseDenoiseFBO;
	VoxelRT::ColorPassFBO ColoredFBO;
	GLClasses::Framebuffer PostProcessingFBO;
	GLClasses::Framebuffer TAAFBO1;
	GLClasses::Framebuffer TAAFBO2;
	GLClasses::TextureArray BlueNoise;
	GLClasses::Framebuffer ShadowFBO_1, ShadowFBO_2;
	GLClasses::Framebuffer ReflectionTraceFBO;
	GLClasses::Framebuffer DownsampledFBO;
	GLClasses::Framebuffer AverageLumaFBO;
	GLClasses::FramebufferRed VolumetricFBO, BlurredVolumetricFBO;
	GLClasses::Framebuffer ReflectionTemporalFBO_1, ReflectionTemporalFBO_2, ReflectionDenoised;

	VoxelRT::BloomFBO BloomFBO(16, 16);

	GLClasses::FramebufferRed SSAOFBO(16, 16);
	GLClasses::FramebufferRed SSAOBlurred(16, 16);

	glm::mat4 CurrentProjection, CurrentView;
	glm::mat4 PreviousProjection, PreviousView;
	glm::mat4 ShadowProjection, ShadowView;
	glm::mat4 ReflectionProjection, ReflectionView;
	glm::vec3 CurrentPosition, PreviousPosition;

	VoxelRT::AtmosphereRenderMap Skymap(64);
	VoxelRT::AtmosphereRenderer AtmosphereRenderer;

	GLClasses::Texture Crosshair;
	GLClasses::Texture BluenoiseTexture;
	Crosshair.CreateTexture("Res/Misc/crosshair.png", false);
	BluenoiseTexture.CreateTexture("Res/Misc/blue_noise.png", false);

	InitialTraceShader.CreateShaderProgramFromFile("Core/Shaders/InitialRayTraceVert.glsl", "Core/Shaders/InitialRayTraceFrag.glsl");
	InitialTraceShader.CompileShaders();
	FinalShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/FBOFrag.glsl");
	FinalShader.CompileShaders();
	DiffuseTraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/DiffuseRayTraceFrag.glsl");
	DiffuseTraceShader.CompileShaders();
	TemporalFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalFilter.glsl");
	TemporalFilter.CompileShaders();
	DenoiseFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SmartDenoise.glsl");
	DenoiseFilter.CompileShaders();
	ColorShader.CreateShaderProgramFromFile("Core/Shaders/ColorPassVert.glsl", "Core/Shaders/ColorPassFrag.glsl");
	ColorShader.CompileShaders();
	PostProcessingShader.CreateShaderProgramFromFile("Core/Shaders/PostProcessingVert.glsl", "Core/Shaders/PostProcessingFrag.glsl");
	PostProcessingShader.CompileShaders();
	TemporalAAShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalAA.glsl");
	TemporalAAShader.CompileShaders();
	ShadowTraceShader.CreateShaderProgramFromFile("Core/Shaders/RayTraceVert.glsl", "Core/Shaders/ShadowRayTraceFrag.glsl");
	ShadowTraceShader.CompileShaders();
	ReflectionTraceShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionTraceFrag.glsl");
	ReflectionTraceShader.CompileShaders();
	SimpleDownsample.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SimpleDownsampleFrag.glsl");
	SimpleDownsample.CompileShaders();
	LumaAverager.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/CalculateAverageLuminance.glsl");
	LumaAverager.CompileShaders();
	VolumetricScattering.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/VolumetricLighting.glsl");
	VolumetricScattering.CompileShaders();
	BilateralBlur.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/BilateralBlur.glsl");
	BilateralBlur.CompileShaders();
	ReflectionDenoiser.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/AdaptiveReflectionDenoise.glsl");
	ReflectionDenoiser.CompileShaders();

	SSAOShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAO.glsl");
	SSAOShader.CompileShaders();
	SSAO_Blur.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAOBlur.glsl");
	SSAO_Blur.CompileShaders();

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

	BloomRenderer::Initialize();
	AverageLumaFBO.SetSize(1, 1);

	///// Set the block texture data uniforms    /////

	InitialTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		std::string name2 = "BLOCK_EMISSIVE_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec4 data;

		float data2 = VoxelRT::BlockDatabase::GetBlockEmissiveTexture(i);
		data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

		InitialTraceShader.SetVector4f(name.c_str(), data);
		InitialTraceShader.SetFloat(name2.c_str(), data2);
	}

	glUseProgram(0);

	ShadowTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec4 data;

		data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

		ShadowTraceShader.SetVector4f(name.c_str(), data);
	}

	glUseProgram(0);

	DiffuseTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec4 data;

		data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.w = VoxelRT::BlockDatabase::GetBlockEmissiveTexture(i);

		DiffuseTraceShader.SetVector4f(name.c_str(), data);
	}

	glUseProgram(0);

	ReflectionTraceShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec4 data;

		data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

		ReflectionTraceShader.SetVector4f(name.c_str(), data);
	}

	glUseProgram(0);

	/////                                     /////

	VoxelRT::InitialRTFBO* InitialTraceFBO = &InitialTraceFBO_1;
	VoxelRT::InitialRTFBO* InitialTraceFBOPrev = &InitialTraceFBO_2;

	glm::vec3 StrongerLightDirection;

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		// Tick the sun and moon
		float time_angle = SunTick * 2.0f;
		glm::mat4 sun_rotation_matrix;

		sun_rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(time_angle), glm::vec3(0.0f, 0.0f, 1.0f));
		SunDirection = glm::vec3(sun_rotation_matrix * glm::vec4(1.0f));
		MoonDirection = glm::vec3(-SunDirection.x, -SunDirection.y, SunDirection.z);
		StrongerLightDirection = -SunDirection.y < 0.01f ? SunDirection : MoonDirection;

		glfwSwapInterval((int)VSync);

		// Resize the framebuffers

		// Diffuse FBOS
		DiffuseTraceFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		DiffuseTemporalFBO1.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		DiffuseTemporalFBO2.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		DiffuseDenoiseFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);

		// MISC
		PostProcessingFBO.SetSize(app.GetWidth(), app.GetHeight());
		ColoredFBO.SetDimensions(app.GetWidth(), app.GetHeight());
		InitialTraceFBO_1.SetDimensions(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution));
		InitialTraceFBO_2.SetDimensions(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution));

		// TAA
		TAAFBO1.SetSize(app.GetWidth(), app.GetHeight());
		TAAFBO2.SetSize(app.GetWidth(), app.GetHeight());

		// 
		DownsampledFBO.SetSize(app.GetWidth() * 0.125f, app.GetHeight() * 0.125f);


		// Reflection and shadow FBOS
		ShadowFBO_1.SetSize(app.GetWidth() * ShadowTraceResolution, app.GetHeight() * ShadowTraceResolution);
		ShadowFBO_2.SetSize(app.GetWidth() * ShadowTraceResolution, app.GetHeight() * ShadowTraceResolution);
		ReflectionTraceFBO.SetSize(app.GetWidth() * ReflectionTraceResolution, app.GetHeight() * ReflectionTraceResolution);
		ReflectionTemporalFBO_1.SetSize(app.GetWidth() * ReflectionTraceResolution, app.GetHeight() * ReflectionTraceResolution);
		ReflectionTemporalFBO_2.SetSize(app.GetWidth() * ReflectionTraceResolution, app.GetHeight() * ReflectionTraceResolution);
		ReflectionDenoised.SetSize(app.GetWidth() * ReflectionTraceResolution, app.GetHeight() * ReflectionTraceResolution);

		// SSAO
		SSAOFBO.SetSize(app.GetWidth() * SSAOResolution, app.GetHeight() * SSAOResolution);
		SSAOBlurred.SetSize(app.GetWidth() * SSAOResolution, app.GetHeight() * SSAOResolution);

		// Volumetrics
		VolumetricFBO.SetSize(app.GetWidth() * VolumetricResolution , app.GetHeight() * VolumetricResolution);
		BlurredVolumetricFBO.SetSize(app.GetWidth() * VolumetricResolution , app.GetHeight() * VolumetricResolution);

		BloomFBO.SetSize(app.GetWidth() * 0.75f, app.GetHeight() * 0.75f);

		///
		GLClasses::Framebuffer& TAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO1 : TAAFBO2;
		GLClasses::Framebuffer& PrevTAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO2 : TAAFBO1;
		GLClasses::Framebuffer& DiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO1 : DiffuseTemporalFBO2;
		GLClasses::Framebuffer& PrevDiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO2 : DiffuseTemporalFBO1;
		GLClasses::Framebuffer& ReflectionTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTemporalFBO_1 : ReflectionTemporalFBO_2;
		GLClasses::Framebuffer& PrevReflectionTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTemporalFBO_2 : ReflectionTemporalFBO_1;
		/// 

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_F2) == GLFW_PRESS)
		{
			InitialTraceShader.Recompile();
			FinalShader.Recompile();
			DiffuseTraceShader.Recompile();
			TemporalFilter.Recompile();
			DenoiseFilter.Recompile();
			ColorShader.Recompile();
			PostProcessingShader.Recompile();
			TemporalAAShader.Recompile();
			ShadowTraceShader.Recompile();
			ReflectionTraceShader.Recompile();
			SSAOShader.Recompile();
			SSAO_Blur.Recompile();
			SimpleDownsample.Recompile();
			LumaAverager.Recompile();
			VolumetricScattering.Recompile();
			BilateralBlur.Recompile();
			ReflectionDenoiser.Recompile();

			BloomRenderer::RecompileShaders();

			///// Set the block texture data uniforms    /////

			InitialTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				std::string name2 = "BLOCK_EMISSIVE_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec4 data;

				float data2 = VoxelRT::BlockDatabase::GetBlockEmissiveTexture(i);
				data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

				InitialTraceShader.SetVector4f(name.c_str(), data);
				InitialTraceShader.SetFloat(name2.c_str(), data2);
			}

			glUseProgram(0);

			ShadowTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec4 data;

				data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

				ShadowTraceShader.SetVector4f(name.c_str(), data);
			}

			glUseProgram(0);

			DiffuseTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec4 data;

				data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.w = VoxelRT::BlockDatabase::GetBlockEmissiveTexture(i);

				DiffuseTraceShader.SetVector4f(name.c_str(), data);
			}

			glUseProgram(0);

			ReflectionTraceShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec4 data;

				data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

				ReflectionTraceShader.SetVector4f(name.c_str(), data);
			}

			glUseProgram(0);

			/////                                     /////

			VoxelRT::Logger::Log("Recompiled!");
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

		if (app.GetCurrentFrame() % 3 == 0)
		{
			AtmosphereRenderer.RenderAtmosphere(Skymap, glm::normalize(SunDirection), 30, 4);
		}

		if (TempView != MainCamera.GetViewMatrix() || ModifiedWorld || app.GetCurrentFrame() % 4 == 0)
		{
			// Swap the initial trace framebuffers
			InitialTraceFBO = InitialTraceFBO == &InitialTraceFBO_1 ? &InitialTraceFBO_2 : &InitialTraceFBO_1;
			InitialTraceFBOPrev = InitialTraceFBO == &InitialTraceFBO_1 ? &InitialTraceFBO_2 : &InitialTraceFBO_1;

			InitialTraceFBO->Bind();
			glDisable(GL_CULL_FACE);
			glDisable(GL_DEPTH_TEST);

			InitialTraceShader.Use();

			InitialTraceShader.SetMatrix4("u_InverseView", inv_view);
			InitialTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
			InitialTraceShader.SetInteger("u_VoxelDataTexture", 0);
			InitialTraceShader.SetInteger("u_AlbedoTextures", 1);
			InitialTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			InitialTraceShader.SetInteger("u_VertCurrentFrame", app.GetCurrentFrame());
			InitialTraceShader.SetVector2f("u_Dimensions", glm::vec2(InitialTraceFBO->GetWidth(), InitialTraceFBO->GetHeight()));
			InitialTraceShader.SetVector2f("u_VertDimensions", glm::vec2(app.GetWidth(), app.GetHeight()));

			/// TEMPORARY ///
			InitialTraceShader.SetInteger("u_GrassBlockProps[0]", VoxelRT::BlockDatabase::GetBlockID("Grass"));
			InitialTraceShader.SetInteger("u_GrassBlockProps[1]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			InitialTraceShader.SetInteger("u_GrassBlockProps[2]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			InitialTraceShader.SetInteger("u_GrassBlockProps[3]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			InitialTraceShader.SetInteger("u_GrassBlockProps[4]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			InitialTraceShader.SetInteger("u_GrassBlockProps[5]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			InitialTraceShader.SetInteger("u_GrassBlockProps[6]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			InitialTraceShader.SetInteger("u_GrassBlockProps[7]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			InitialTraceShader.SetInteger("u_GrassBlockProps[8]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			InitialTraceShader.SetInteger("u_GrassBlockProps[9]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			InitialTraceFBO->Unbind();
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
		DiffuseTraceShader.SetInteger("u_BlockPBRTextures", 8);
		DiffuseTraceShader.SetInteger("u_BlockEmissiveTextures", 11);
		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
		DiffuseTraceShader.SetMatrix4("u_InverseView", inv_view);
		DiffuseTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
		DiffuseTraceShader.SetVector2f("u_Dimensions", glm::vec2(DiffuseTraceFBO.GetWidth(), DiffuseTraceFBO.GetHeight()));
		DiffuseTraceShader.SetFloat("u_Time", glfwGetTime() * 1.2f);
		DiffuseTraceShader.SetFloat("u_DiffuseLightIntensity", DiffuseLightIntensity);
		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
		DiffuseTraceShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
		DiffuseTraceShader.SetVector3f("u_SunDirection", SunDirection);
		DiffuseTraceShader.SetVector3f("u_MoonDirection", MoonDirection);

		DiffuseTraceShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
		DiffuseTraceShader.SetMatrix4("u_ShadowView", ShadowView);
		DiffuseTraceShader.SetInteger("u_ShadowMap", 9);
		DiffuseTraceShader.SetInteger("u_BlueNoiseTexture", 10);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetTexture());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetDataTexture());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlueNoise.GetTextureArray());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, ShadowFBO_1.GetTexture());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTraceFBO.Unbind();

		///// Temporal filter the calculated diffuse /////

		DiffuseTemporalFBO.Bind();
		TemporalFilter.Use();

		TemporalFilter.SetInteger("u_CurrentColorTexture", 0);
		TemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
		TemporalFilter.SetInteger("u_PreviousColorTexture", 2);
		TemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);

		TemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
		TemporalFilter.SetMatrix4("u_View", CurrentView);
		TemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
		TemporalFilter.SetMatrix4("u_PrevView", PreviousView);

		TemporalFilter.SetFloat("u_MixModifier", 0.8f);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTraceFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetPositionTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTemporalFBO.Unbind();

		//// DENOISE TEMPORALLY FILTERED OUTPUT ////

		// ---- Denoise the diffuse ----

		DiffuseDenoiseFBO.Bind();
		DenoiseFilter.Use();

		DenoiseFilter.SetInteger("u_Texture", 0);
		DenoiseFilter.SetInteger("u_NormalTexture", 1);
		DenoiseFilter.SetInteger("u_InitialTracePositionTexture", 2);
		DenoiseFilter.SetInteger("u_Radius", 12);
		DenoiseFilter.SetFloat("u_EdgeThreshold", 0.5f);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseDenoiseFBO.Unbind();


		// ---- SHADOW TRACE ----
		GLClasses::Framebuffer& ShadowFBO = app.GetCurrentFrame() % 2 == 0 ? ShadowFBO_1 : ShadowFBO_2;
		GLClasses::Framebuffer& PrevShadowFBO = app.GetCurrentFrame() % 2 == 0 ? ShadowFBO_2 : ShadowFBO_1;

		{
			bool DoFullTrace = FullyDynamicShadows ? true : (((app.GetCurrentFrame() % 5 == 0) || ModifiedWorld) ? true : false);

			ShadowFBO.Bind();
			ShadowTraceShader.Use();

			ShadowTraceShader.SetInteger("u_PositionTexture", 0);
			ShadowTraceShader.SetInteger("u_VoxelData", 1);
			ShadowTraceShader.SetInteger("u_AlbedoTextures", 2);
			ShadowTraceShader.SetInteger("u_NormalTexture", 3);
			ShadowTraceShader.SetInteger("u_PrevShadowFBO", 4);
			ShadowTraceShader.SetVector3f("u_LightDirection", StrongerLightDirection);
			ShadowTraceShader.SetVector3f("u_PlayerPosition", MainCamera.GetPosition());
			ShadowTraceShader.SetBool("u_DoFullTrace", DoFullTrace);
			ShadowTraceShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
			ShadowTraceShader.SetMatrix4("u_ShadowView", ShadowView);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D, PrevShadowFBO.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ShadowFBO.Unbind();

			ShadowProjection = CurrentProjection;
			ShadowView = CurrentView;
		}

		// ---- REFLECTION TRACE ----


		{
			ReflectionTraceFBO.Bind();
			ReflectionTraceShader.Use();

			ReflectionTraceShader.SetInteger("u_PositionTexture", 0);
			ReflectionTraceShader.SetInteger("u_PBRTexture", 1);
			ReflectionTraceShader.SetInteger("u_InitialTraceNormalTexture", 2);
			ReflectionTraceShader.SetInteger("u_BlockNormalTextures", 4);
			ReflectionTraceShader.SetInteger("u_BlockAlbedoTextures", 5);
			ReflectionTraceShader.SetInteger("u_BlockPBRTextures", 6);
			ReflectionTraceShader.SetInteger("u_Skymap", 7);
			ReflectionTraceShader.SetInteger("u_VoxelData", 8);
			ReflectionTraceShader.SetInteger("u_BlueNoiseTexture", 9);
			ReflectionTraceShader.SetFloat("u_ReflectionTraceRes", ReflectionTraceResolution);
			ReflectionTraceShader.SetVector3f("u_SunDirection", SunDirection);
			ReflectionTraceShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
			ReflectionTraceShader.SetVector2f("u_Dimensions", glm::vec2(ReflectionTraceFBO.GetWidth(), ReflectionTraceFBO.GetHeight()));
			ReflectionTraceShader.SetFloat("u_Time", glfwGetTime());
			ReflectionTraceShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());

			ReflectionTraceShader.SetInteger("u_GrassBlockProps[0]", VoxelRT::BlockDatabase::GetBlockID("Grass"));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[1]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[2]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[3]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[4]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[5]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[6]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[7]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[8]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			ReflectionTraceShader.SetInteger("u_GrassBlockProps[9]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));


			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetDataTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE6);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

			glActiveTexture(GL_TEXTURE7);
			glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetTexture());

			glActiveTexture(GL_TEXTURE8);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE9);
			glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionTraceFBO.Unbind();
		}

		// Temporally filter it
		{

			ReflectionTemporalFBO.Bind();
			TemporalFilter.Use();

			TemporalFilter.SetInteger("u_CurrentColorTexture", 0);
			TemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
			TemporalFilter.SetInteger("u_PreviousColorTexture", 2);
			TemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);

			TemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
			TemporalFilter.SetMatrix4("u_View", CurrentView);
			TemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
			TemporalFilter.SetMatrix4("u_PrevView", PreviousView);

			TemporalFilter.SetFloat("u_MixModifier", ModifiedWorld ? 0.5f : 0.8f); // Brutal temporal filtering

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, PrevReflectionTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetPositionTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionTemporalFBO.Unbind();
		}

		// Denoise reflection trace
		{
			ReflectionDenoised.Bind();
			ReflectionDenoiser.Use();

			ReflectionDenoiser.SetInteger("u_Texture", 0);
			ReflectionDenoiser.SetInteger("u_PBRTexture", 1);
			ReflectionDenoiser.SetInteger("u_Radius", 12);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetPBRTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionDenoised.Unbind();
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
		ColorShader.SetInteger("u_BlockEmissiveTextures", 11);
		ColorShader.SetMatrix4("u_InverseView", inv_view);
		ColorShader.SetMatrix4("u_InverseProjection", inv_projection);
		ColorShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
		ColorShader.SetMatrix4("u_ShadowView", ShadowView);
		ColorShader.SetMatrix4("u_ReflectionProjection", ReflectionProjection);
		ColorShader.SetMatrix4("u_ReflectionView", ReflectionView);
		ColorShader.SetVector2f("u_InitialTraceResolution", glm::vec2(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution)));
		ColorShader.SetVector3f("u_SunDirection", SunDirection);
		ColorShader.SetVector3f("u_MoonDirection", MoonDirection);
		ColorShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
		ColorShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
		ColorShader.SetFloat("u_Time", glfwGetTime());

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetDataTexture());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetTexture());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, ShadowFBO.GetTexture());

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlueNoise.GetTextureArray());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, ReflectionDenoised.GetTexture());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		ColoredFBO.Unbind();

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
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevTAAFBO.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		// ---- SSAO ---- //

		if (SSAO)
		{
			glDisable(GL_DEPTH_TEST);
			glDisable(GL_CULL_FACE);
			glDisable(GL_BLEND);

			SSAOShader.Use();
			SSAOFBO.Bind();

			SSAOShader.SetInteger("u_PositionTexture", 0);
			SSAOShader.SetInteger("u_NormalTexture", 1);
			SSAOShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			SSAOShader.SetInteger("SAMPLE_SIZE", 20);
			SSAOShader.SetMatrix4("u_ViewMatrix", MainCamera.GetViewMatrix());
			SSAOShader.SetMatrix4("u_ProjectionMatrix", MainCamera.GetProjectionMatrix());
			SSAOShader.SetVector2f("u_Dimensions", glm::vec2(SSAOFBO.GetWidth(), SSAOFBO.GetHeight()));
			SSAOShader.SetFloat("u_Time", glfwGetTime());

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetNormalTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
			SSAOFBO.Unbind();
			glUseProgram(0);

			SSAOBlurred.Bind();
			SSAO_Blur.Use();
			SSAO_Blur.SetInteger("u_Texture", 0);
			SSAO_Blur.SetVector2f("u_SketchSize", glm::vec2(SSAOFBO.GetWidth(), SSAOFBO.GetHeight()));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, SSAOFBO.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			SSAOBlurred.Unbind();

			glUseProgram(0);
		}

		// ---- Bloom ----
		
		if (Bloom)
		{
			BloomRenderer::RenderBloom(BloomFBO, ColoredFBO.GetColorTexture());
		}

		// ---- Auto Exposure ----

		float ComputedExposure = 3.0f;

		if (AutoExposure)
		{
			DownsampledFBO.Bind();
			SimpleDownsample.Use();
			SimpleDownsample.SetInteger("u_Texture", 0);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, TAAFBO.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			AverageLumaFBO.Bind();
			LumaAverager.Use();

			LumaAverager.SetInteger("u_Texture", 0);
			
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, DownsampledFBO.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			AverageLumaFBO.Bind();
			GLfloat avg_col[3];
			glReadPixels(0, 0, 1, 1, GL_RGB, GL_FLOAT, &avg_col);
			GLfloat lum = avg_col[0];
			ComputedExposure = lum * 10.0f;
			ComputedExposure = glm::clamp(ComputedExposure, 0.5f, 3.6f);

			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}

		// ---- Volumetric Scattering ----

		if (GodRays)
		{
			VolumetricFBO.Bind();
			VolumetricScattering.Use();

			VolumetricScattering.SetInteger("u_PositionTexture", 0);
			VolumetricScattering.SetInteger("u_BlueNoiseTexture", 1);
			VolumetricScattering.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
			VolumetricScattering.SetVector2f("u_Dimensions", glm::vec2(VolumetricFBO.GetWidth(), VolumetricFBO.GetHeight()));
			VolumetricScattering.SetMatrix4("u_ProjectionMatrix", MainCamera.GetProjectionMatrix());
			VolumetricScattering.SetMatrix4("u_ViewMatrix", MainCamera.GetViewMatrix());

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			VolumetricFBO.Unbind();

			// Horizontal Blur 

			BlurredVolumetricFBO.Bind();
			BilateralBlur.Use();

			BilateralBlur.SetFloat("u_Sharpness", 0.5f);
			BilateralBlur.SetVector2f("u_InvResolutionDirection", glm::vec2(1.0f / BlurredVolumetricFBO.GetWidth(), 0.0f));
			BilateralBlur.SetInteger("u_ColorTexture", 0);
			BilateralBlur.SetInteger("u_PositionTexture", 1);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, VolumetricFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			BlurredVolumetricFBO.Unbind();

			// Vertical blur

			VolumetricFBO.Bind();

			BilateralBlur.Use();

			BilateralBlur.SetFloat("u_Sharpness", 0.5f);
			BilateralBlur.SetVector2f("u_InvResolutionDirection", glm::vec2(0.0f, 1.0f / BlurredVolumetricFBO.GetHeight()));
			BilateralBlur.SetInteger("u_ColorTexture", 0);
			BilateralBlur.SetInteger("u_PositionTexture", 1);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, BlurredVolumetricFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			VolumetricFBO.Unbind();
		}


		// ---- POST PROCESSING ----

		PostProcessingFBO.Bind();
		PostProcessingShader.Use();

		PostProcessingShader.SetMatrix4("u_InverseView", inv_view);
		PostProcessingShader.SetMatrix4("u_InverseProjection", inv_projection);
		PostProcessingShader.SetInteger("u_FramebufferTexture", 0);
		PostProcessingShader.SetInteger("u_PositionTexture", 1);
		PostProcessingShader.SetInteger("u_BlueNoise", 2);
		PostProcessingShader.SetInteger("u_SSAOTexture", 3);
		PostProcessingShader.SetInteger("u_PositionTexture", 4);
		PostProcessingShader.SetInteger("u_VolumetricTexture", 10);
		PostProcessingShader.SetInteger("u_GodRaysStepCount", GodRaysStepCount);
		PostProcessingShader.SetVector3f("u_SunDirection", SunDirection);
		PostProcessingShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
		PostProcessingShader.SetVector2f("u_Dimensions", glm::vec2(PostProcessingFBO.GetWidth(), PostProcessingFBO.GetHeight()));
		PostProcessingShader.SetMatrix4("u_ProjectionMatrix", MainCamera.GetProjectionMatrix());
		PostProcessingShader.SetMatrix4("u_ViewMatrix", MainCamera.GetViewMatrix());
		PostProcessingShader.SetBool("u_SunIsStronger", StrongerLightDirection == SunDirection);
		PostProcessingShader.SetBool("u_LensFlare", LensFlare);
		PostProcessingShader.SetBool("u_GodRays", GodRays);
		PostProcessingShader.SetBool("u_SSAO", SSAO);
		PostProcessingShader.SetBool("u_Bloom", Bloom);
		PostProcessingShader.SetBool("u_SSGodRays", FakeGodRays);
		PostProcessingShader.SetFloat("u_LensFlareIntensity", LensFlareIntensity);
		PostProcessingShader.SetFloat("u_Exposure", ComputedExposure);

		// Set the bloom mips
		PostProcessingShader.SetInteger("u_BloomMips[0]", 5);
		PostProcessingShader.SetInteger("u_BloomMips[1]", 6);
		PostProcessingShader.SetInteger("u_BloomMips[2]", 7);
		PostProcessingShader.SetInteger("u_BloomMips[3]", 8);
		PostProcessingShader.SetInteger("u_ShadowTexture", 9);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TAAFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, SSAOBlurred.GetTexture());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		// Bloom mips
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mip0);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mip1);

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mip2);

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mip3);
		//

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, ShadowFBO.GetTexture());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, VolumetricFBO.GetTexture());

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
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		RendererUI.RenderQuad(glm::vec2(floor((float)app.GetWidth() / 2.0f), floor((float)app.GetHeight() / 2.0f)), &Crosshair, &OCamera);

		// Finish Frame
		app.FinishFrame();
		
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glUseProgram(0);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		std::string title = "Voxel RT | "; title += BlockDatabase::GetBlockName(world->GetCurrentBlock()); title += "  ";
		GLClasses::DisplayFrameRate(app.GetWindow(), title);

		ModifiedWorld = false;
	}

	SaveWorld(world, world_name);
	delete world;

	return;
}
