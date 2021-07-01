#include "Pipeline.h"


#include <chrono>



static VoxelRT::Player MainPlayer;
static bool VSync = false;

static bool CloudsEnabled = true;
static float CloudCoverage = 0.135f;
static bool CloudBayer = true;
static float CloudDetailContribution = 0.01f;
static bool CloudHighQuality = false;

static bool DoSecondSpatialPass = false;

static float InitialTraceResolution = 0.500f;
static float DiffuseTraceResolution = 0.200f; // 1/5th res + 4 spp = 0.8 spp

static float ShadowTraceResolution = 0.750;
static float ReflectionTraceResolution = 0.300; // Resolution is key for specular indirect.
static float SSAOResolution = 0.35f;
static float RTAOResolution = 0.125f;

static float VolumetricResolution = 0.5f;

static float SunTick = 50.0f;
static float DiffuseLightIntensity = 80.0f;
static float LensFlareIntensity = 0.075f;
static float BloomQuality = 1.0f;

static int DiffuseSPP = 4;
static int ReflectionSPP = 4;

static bool TAA = true;
static bool Bloom = true;

static bool GodRays = false;
static bool FakeGodRays = false;
static bool RoughReflections = true;
static bool DenoiseReflections = true;
static bool RenderParticles = true;

static bool LensFlare = false;
static bool SSAO = false;
static bool RTAO = false;
static bool POM = false;
static bool HighQualityPOM = false;

static bool CheckerboardClouds = true;

static bool FullyDynamicShadows = true;

static int GodRaysStepCount = 12;



static bool AutoExposure = false;
static bool ExponentialFog = false;

static glm::vec3 SunDirection;
static glm::vec3 MoonDirection;

static VoxelRT::World* world = nullptr;
static bool ModifiedWorld = false;
static VoxelRT::FPSCamera& MainCamera = MainPlayer.Camera;
static VoxelRT::OrthographicCamera OCamera(0.0f, 800.0f, 0.0f, 600.0f);

static float Frametime;
static float DeltaTime;

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
		if (ImGui::Begin("Settings"))
		{
			ImGui::Text("Player Position : %f, %f, %f", MainCamera.GetPosition().x, MainCamera.GetPosition().y, MainCamera.GetPosition().z);
			ImGui::Text("Camera Front : %f, %f, %f", MainCamera.GetFront().x, MainCamera.GetFront().y, MainCamera.GetFront().z);
			ImGui::SliderFloat("Initial Trace Resolution", &InitialTraceResolution, 0.1f, 1.0f);
			ImGui::SliderFloat("Diffuse Trace Resolution ", &DiffuseTraceResolution, 0.1f, 1.25f);
			ImGui::SliderFloat("Shadow Trace Resolution ", &ShadowTraceResolution, 0.1f, 1.25f);
			ImGui::SliderFloat("Reflection Trace Resolution ", &ReflectionTraceResolution, 0.1f, 0.8f);
			ImGui::SliderFloat("SSAO Render Resolution ", &SSAOResolution, 0.1f, 0.9f);
			ImGui::SliderFloat("Diffuse Light Intensity ", &DiffuseLightIntensity, -10.0, 80.0f);
			ImGui::SliderFloat("Sun Time ", &SunTick, 0.1f, 256.0f);
			ImGui::SliderFloat("Lens Flare Intensity ", &LensFlareIntensity, 0.05f, 1.25f);
			ImGui::SliderFloat("RTAO Resolution ", &RTAOResolution, 0.1f, 0.9f);
			ImGui::SliderFloat("Bloom Quality ", &BloomQuality, 0.1f, 1.5f);
			ImGui::SliderInt("God ray raymarch step count", &GodRaysStepCount, 8, 64);
			ImGui::SliderInt("Diffuse Trace SPP", &DiffuseSPP, 1, 32);
			ImGui::SliderInt("Reflection Trace SPP", &ReflectionSPP, 1, 64);
			ImGui::Checkbox("Do second spatial filtering pass (For indirect, more expensive, reduces noise) ?", &DoSecondSpatialPass);
			ImGui::Checkbox("Rough reflections?", &RoughReflections);
			ImGui::Checkbox("Denoise reflections?", &DenoiseReflections);
			ImGui::Checkbox("Particles?", &RenderParticles);
			ImGui::Checkbox("Fully Dynamic Shadows? (Fixes shadow artifacts)", &FullyDynamicShadows);
			ImGui::Checkbox("Ray traced ambient occlusion (Slower, more accurate)?", &RTAO);
			ImGui::Checkbox("POM? (WIP)", &POM);
			ImGui::Checkbox("High Quality POM?", &HighQualityPOM);
			ImGui::Checkbox("Temporal Anti Aliasing", &TAA);
			ImGui::Checkbox("Volumetric Clouds?", &CloudsEnabled);
			ImGui::Checkbox("High Quality Clouds? (Doubles the ray march step count)", &CloudHighQuality);
			ImGui::Checkbox("Use Bayer Dither for clouds? (Uses white noise if disabled)", &CloudBayer);
			ImGui::SliderFloat("Volumetric Cloud Coverage", &CloudCoverage, 0.01f, 0.6f);
			ImGui::SliderFloat("Volumetric Cloud Detail Contribution", &CloudDetailContribution, 0.0f, 1.5f);

			ImGui::Checkbox("Checkerboard clouds?", &CheckerboardClouds);
			ImGui::Checkbox("Lens Flare?", &LensFlare);
			ImGui::Checkbox("(Implementation - 1) God Rays? (Slower)", &GodRays);
			ImGui::Checkbox("(Implementation - 2) God Rays? (faster, more crisp, Adjust the step count in the menu)", &FakeGodRays);
			ImGui::Checkbox("Screen Space Ambient Occlusion?", &SSAO);
			ImGui::Checkbox("Exponential Fog?", &ExponentialFog);
			ImGui::Checkbox("Bloom (Expensive!) ?", &Bloom);
			ImGui::Checkbox("Auto Exposure (Very very WIP!) ?", &AutoExposure);
		} ImGui::End();

		if (ImGui::Begin("Other Settings"))
		{
			ImGui::SliderFloat("Mouse Sensitivity", &MainPlayer.Sensitivity, 0.025f, 1.0f);
			ImGui::SliderFloat("Player Speed", &MainPlayer.Speed, 0.025f, 1.0f);
			ImGui::Checkbox("VSync", &VSync);

			if (ImGui::Button("Reset"))
			{
				MainPlayer.Sensitivity = 0.25f;
				MainPlayer.Speed = 0.045f;
				VSync = false;
			}
		} ImGui::End();
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
				world->Raycast(0, MainCamera.GetPosition(), MainCamera.GetFront());
				ModifiedWorld = true;
			}
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_RIGHT && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(1, MainCamera.GetPosition(), MainCamera.GetFront());
				ModifiedWorld = true;
			}
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_MIDDLE && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(2, MainCamera.GetPosition(), MainCamera.GetFront());
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
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

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
	world->InitializeDistanceGenerator();
	world->GenerateDistanceField();

	VoxelRT::Renderer2D RendererUI;

	GLClasses::VertexBuffer VBO;
	GLClasses::VertexArray VAO;
	GLClasses::Shader InitialTraceShader;
	GLClasses::Shader FinalShader;
	GLClasses::Shader DiffuseTraceShader;
	GLClasses::Shader MainTemporalFilter;
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
	GLClasses::Shader RTAOShader;
	GLClasses::Shader SpatialFilter;
	GLClasses::Shader SpatialCleanup;

	VoxelRT::InitialRTFBO InitialTraceFBO_1;
	VoxelRT::InitialRTFBO InitialTraceFBO_2;
	GLClasses::Framebuffer DiffuseTraceFBO;
	GLClasses::Framebuffer DiffuseTemporalFBO1;
	GLClasses::Framebuffer DiffuseTemporalFBO2;
	GLClasses::Framebuffer DiffuseDenoiseFBO, DiffuseDenoisedFBO2;
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
	GLClasses::Framebuffer RTAO_FBO, RTAO_TemporalFBO_1, RTAO_TemporalFBO_2;

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
	MainTemporalFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalFilter.glsl");
	MainTemporalFilter.CompileShaders();
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
	RTAOShader.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/RaytracedAO.glsl");
	RTAOShader.CompileShaders();
	SpatialFilter.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SpatialFilter.glsl");
	SpatialFilter.CompileShaders();
	SpatialCleanup.CreateShaderProgramFromFile("Core/Shaders/FBOVert.glsl", "Core/Shaders/SpatialFinalCleanup.glsl");
	SpatialCleanup.CompileShaders();

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

	RTAOShader.Use();

	for (int i = 0; i < 128; i++)
	{
		// BLOCK_TEXTURE_DATA

		std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
		glm::vec4 data;

		data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
		data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

		RTAOShader.SetVector4f(name.c_str(), data);
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
		data.w = (float)VoxelRT::BlockDatabase::GetBlockEmissiveTexture(i);

		ReflectionTraceShader.SetVector4f(name.c_str(), data);
	}

	glUseProgram(0);

	Clouds::CloudRenderer::Initialize();

	/////                                     /////

	VoxelRT::InitialRTFBO* InitialTraceFBO = &InitialTraceFBO_1;
	VoxelRT::InitialRTFBO* InitialTraceFBOPrev = &InitialTraceFBO_2;

	glm::vec3 StrongerLightDirection;


	GLfloat PreviousLuma = 3.0f;

	float CameraExposure = 1.0f;
	float PrevCameraExposure = 1.0f;

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
		float DiffuseTraceResolution2 = DiffuseTraceResolution + 0.125f;
		DiffuseTraceFBO.SetSize(app.GetWidth() * DiffuseTraceResolution, app.GetHeight() * DiffuseTraceResolution);
		DiffuseTemporalFBO1.SetSize(app.GetWidth() * DiffuseTraceResolution2, app.GetHeight() * DiffuseTraceResolution2);
		DiffuseTemporalFBO2.SetSize(app.GetWidth() * DiffuseTraceResolution2, app.GetHeight() * DiffuseTraceResolution2);
		DiffuseDenoiseFBO.SetSize(app.GetWidth() * DiffuseTraceResolution2, app.GetHeight() * DiffuseTraceResolution2);
		DiffuseDenoisedFBO2.SetSize(app.GetWidth() * DiffuseTraceResolution2, app.GetHeight() * DiffuseTraceResolution2);

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

		BloomFBO.SetSize(app.GetWidth() * BloomQuality, app.GetHeight() * BloomQuality);

		//

		float RTAO_Res2 = glm::max(RTAOResolution, 0.5f);
		RTAO_FBO.SetSize(app.GetWidth() * RTAOResolution, app.GetHeight() * RTAOResolution);
		RTAO_TemporalFBO_1.SetSize(app.GetWidth() * RTAO_Res2, app.GetHeight() * RTAO_Res2);
		RTAO_TemporalFBO_2.SetSize(app.GetWidth() * RTAO_Res2, app.GetHeight() * RTAO_Res2);

		///
		GLClasses::Framebuffer& TAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO1 : TAAFBO2;
		GLClasses::Framebuffer& PrevTAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO2 : TAAFBO1;
		GLClasses::Framebuffer& DiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO1 : DiffuseTemporalFBO2;
		GLClasses::Framebuffer& PrevDiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO2 : DiffuseTemporalFBO1;
		GLClasses::Framebuffer& ReflectionTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTemporalFBO_1 : ReflectionTemporalFBO_2;
		GLClasses::Framebuffer& PrevReflectionTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTemporalFBO_2 : ReflectionTemporalFBO_1;
		GLClasses::Framebuffer& RTAOTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? RTAO_TemporalFBO_1: RTAO_TemporalFBO_2;
		GLClasses::Framebuffer& PrevRTAOTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? RTAO_TemporalFBO_2 : RTAO_TemporalFBO_1;
		/// 

		if (glfwGetKey(app.GetWindow(), GLFW_KEY_F2) == GLFW_PRESS)
		{
			InitialTraceShader.Recompile();
			FinalShader.Recompile();
			DiffuseTraceShader.Recompile();
			MainTemporalFilter.Recompile();
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
			RTAOShader.Recompile();
			SpatialFilter.Recompile();
			SpatialCleanup.Recompile();

			world->m_ParticleEmitter.Recompile();
			Clouds::CloudRenderer::RecompileShaders();
			BloomRenderer::RecompileShaders();
			AtmosphereRenderer.Recompile();

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

			RTAOShader.Use();

			for (int i = 0; i < 128; i++)
			{
				// BLOCK_TEXTURE_DATA

				std::string name = "BLOCK_TEXTURE_DATA[" + std::to_string(i) + "]";
				glm::vec4 data;

				data.x = VoxelRT::BlockDatabase::GetBlockTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.y = VoxelRT::BlockDatabase::GetBlockNormalTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.z = VoxelRT::BlockDatabase::GetBlockPBRTexture(i, VoxelRT::BlockDatabase::BlockFaceType::Top);
				data.w = VoxelRT::BlockDatabase::IsBlockTransparent(i);

				RTAOShader.SetVector4f(name.c_str(), data);
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
				data.w = VoxelRT::BlockDatabase::GetBlockEmissiveTexture(i);

				ReflectionTraceShader.SetVector4f(name.c_str(), data);
			}

			glUseProgram(0);

			/////                                     /////

			VoxelRT::Logger::Log("Recompiled!");
		}

		MainPlayer.OnUpdate(app.GetWindow(), world, DeltaTime * 6.9f);
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
			InitialTraceShader.SetInteger("u_DistanceFieldTexture", 2);
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

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

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
		DiffuseTraceShader.SetInteger("u_DistanceFieldTexture", 13);

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
		DiffuseTraceShader.SetInteger("u_SPP", DiffuseSPP);

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

		glActiveTexture(GL_TEXTURE13);
		glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseTraceFBO.Unbind();

		///// Temporal filter the calculated diffuse /////

		DiffuseTemporalFBO.Bind();
		MainTemporalFilter.Use();

		MainTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
		MainTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
		MainTemporalFilter.SetInteger("u_PreviousColorTexture", 2);
		MainTemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);

		MainTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
		MainTemporalFilter.SetMatrix4("u_View", CurrentView);
		MainTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
		MainTemporalFilter.SetMatrix4("u_PrevView", PreviousView);

		MainTemporalFilter.SetFloat("u_MinimumMix", 0.3f);
		MainTemporalFilter.SetFloat("u_MaximumMix", 0.975f);
		MainTemporalFilter.SetInteger("u_TemporalQuality", 1);
		MainTemporalFilter.SetBool("u_ReflectionTemporal", false);

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

		// Spatial pass 1
		{
			DiffuseDenoisedFBO2.Bind();
			SpatialFilter.Use();

			SpatialFilter.SetInteger("u_InputTexture", 0);
			SpatialFilter.SetInteger("u_PositionTexture", 1);
			SpatialFilter.SetInteger("u_NormalTexture", 2);
			SpatialFilter.SetInteger("u_Step", 1);
			SpatialFilter.SetBool("u_Dir", true);
			SpatialFilter.SetVector2f("u_Dimensions", glm::vec2(DiffuseTemporalFBO.GetWidth(), DiffuseTemporalFBO.GetHeight()));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}

		{
			DiffuseDenoiseFBO.Bind();
			SpatialFilter.Use();

			SpatialFilter.SetInteger("u_InputTexture", 0);
			SpatialFilter.SetInteger("u_PositionTexture", 1);
			SpatialFilter.SetInteger("u_NormalTexture", 2);
			SpatialFilter.SetInteger("u_Step", 1);
			SpatialFilter.SetBool("u_Dir", false);
			SpatialFilter.SetVector2f("u_Dimensions", glm::vec2(DiffuseDenoisedFBO2.GetWidth(), DiffuseDenoisedFBO2.GetHeight()));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, DiffuseDenoisedFBO2.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}

		// Spatial pass.. 2 ?

		if (DoSecondSpatialPass)
		{
			{
				DiffuseDenoisedFBO2.Bind();
				SpatialFilter.Use();

				SpatialFilter.SetInteger("u_InputTexture", 0);
				SpatialFilter.SetInteger("u_PositionTexture", 1);
				SpatialFilter.SetInteger("u_NormalTexture", 2);
				SpatialFilter.SetInteger("u_Step", 1);
				SpatialFilter.SetBool("u_Dir", true);
				SpatialFilter.SetVector2f("u_Dimensions", glm::vec2(DiffuseDenoiseFBO.GetWidth(), DiffuseDenoiseFBO.GetHeight()));

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}

			{
				DiffuseDenoiseFBO.Bind();
				SpatialFilter.Use();

				SpatialFilter.SetInteger("u_InputTexture", 0);
				SpatialFilter.SetInteger("u_PositionTexture", 1);
				SpatialFilter.SetInteger("u_NormalTexture", 2);
				SpatialFilter.SetInteger("u_Step", 1);
				SpatialFilter.SetBool("u_Dir", false);
				SpatialFilter.SetVector2f("u_Dimensions", glm::vec2(DiffuseDenoisedFBO2.GetWidth(), DiffuseDenoisedFBO2.GetHeight()));

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, DiffuseDenoisedFBO2.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}
		}

		glUseProgram(0); 
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

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
			ShadowTraceShader.SetInteger("u_DistanceFieldTexture", 5);

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

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

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
			ReflectionTraceShader.SetInteger("u_DistanceFieldTexture", 10);
			ReflectionTraceShader.SetInteger("u_BlockEmissiveTextures", 11);
			ReflectionTraceShader.SetFloat("u_ReflectionTraceRes", ReflectionTraceResolution);
			ReflectionTraceShader.SetVector3f("u_SunDirection", SunDirection);
			ReflectionTraceShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
			ReflectionTraceShader.SetVector2f("u_Dimensions", glm::vec2(ReflectionTraceFBO.GetWidth(), ReflectionTraceFBO.GetHeight()));
			ReflectionTraceShader.SetFloat("u_Time", glfwGetTime());
			ReflectionTraceShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
			ReflectionTraceShader.SetBool("u_RoughReflections", RoughReflections);
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
			ReflectionTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			
			ReflectionTraceShader.SetInteger("u_SPP", ReflectionSPP);

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

			glActiveTexture(GL_TEXTURE10);
			glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE11);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionTraceFBO.Unbind();
		}

		// Temporally filter it
		if (RoughReflections)
		{
			ReflectionTemporalFBO.Bind();
			MainTemporalFilter.Use();

			MainTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
			MainTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
			MainTemporalFilter.SetInteger("u_PreviousColorTexture", 2);
			MainTemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);

			MainTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
			MainTemporalFilter.SetMatrix4("u_View", CurrentView);
			MainTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
			MainTemporalFilter.SetMatrix4("u_PrevView", PreviousView);

			MainTemporalFilter.SetFloat("u_MinimumMix", 0.0f); // Brutal temporal filtering
			MainTemporalFilter.SetFloat("u_MaximumMix", 0.95f);
			MainTemporalFilter.SetInteger("u_TemporalQuality", 1);
			MainTemporalFilter.SetBool("u_ReflectionTemporal", true);

			MainTemporalFilter.SetVector3f("u_PrevCameraPos", PreviousPosition);
			MainTemporalFilter.SetVector3f("u_CurrentCameraPos", MainCamera.GetPosition());

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
		if (DenoiseReflections && RoughReflections)
		{
			ReflectionDenoised.Bind();
			ReflectionDenoiser.Use();
		
			ReflectionDenoiser.SetInteger("u_Texture", 0);
			ReflectionDenoiser.SetInteger("u_PBRTexture", 1);
			ReflectionDenoiser.SetInteger("u_Radius", 5);
			ReflectionDenoiser.SetBool("u_Moved", PreviousPosition != MainCamera.GetPosition());
		
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture());
		
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetPBRTexture());
		
			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		
			ReflectionDenoised.Unbind();
		}

		// ---- RT AO ----

		if (RTAO)
		{
			RTAO_FBO.Bind();
			RTAOShader.Use();

			RTAOShader.SetInteger("u_PositionTexture", 0);
			RTAOShader.SetInteger("u_VoxelData", 1);
			RTAOShader.SetInteger("u_NormalTexture", 2);
			RTAOShader.SetInteger("u_BlockAlbedoTextures", 3);
			RTAOShader.SetInteger("u_BlockNormalTextures", 4);
			RTAOShader.SetInteger("u_DataTexture", 5);
			RTAOShader.SetFloat("u_Time", glfwGetTime());

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetNormalTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetDataTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			RTAO_FBO.Unbind();

			// Temporally filter it
			RTAOTemporalFBO.Bind();
			MainTemporalFilter.Use();

			MainTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
			MainTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
			MainTemporalFilter.SetInteger("u_PreviousColorTexture", 2);
			MainTemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);

			MainTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
			MainTemporalFilter.SetMatrix4("u_View", CurrentView);
			MainTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
			MainTemporalFilter.SetMatrix4("u_PrevView", PreviousView);

			MainTemporalFilter.SetFloat("u_MinimumMix", 0.25f); 
			MainTemporalFilter.SetFloat("u_MaximumMix", 0.975f); 
			MainTemporalFilter.SetInteger("u_TemporalQuality", 1); 
			MainTemporalFilter.SetBool("u_ReflectionTemporal", false);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, RTAO_FBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, PrevRTAOTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetPositionTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			RTAOTemporalFBO.Unbind();
		}

		// ---- RENDER CLOUDS ----
		GLuint CloudData = 0;
		if (CloudsEnabled)
		{
			CloudData = Clouds::CloudRenderer::Update(MainCamera, PreviousProjection,
				PreviousView, CurrentPosition,
				PreviousPosition, VAO, StrongerLightDirection, BluenoiseTexture.GetTextureID(),
				app.GetWidth(), app.GetHeight(), app.GetCurrentFrame(), Skymap.GetTexture(), PreviousPosition);

			Clouds::CloudRenderer::SetChecker(CheckerboardClouds);
			Clouds::CloudRenderer::SetCoverage(CloudCoverage);
			Clouds::CloudRenderer::SetBayer(CloudBayer);
			Clouds::CloudRenderer::SetDetailContribution(CloudDetailContribution);
			Clouds::CloudRenderer::SetQuality(CloudHighQuality);
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
		ColorShader.SetInteger("u_CloudData", 12);
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
		ColorShader.SetFloat("u_GrassblockAlbedoID", BlockDatabase::GetBlockTexture("Grass", BlockDatabase::BlockFaceType::Front));
		ColorShader.SetFloat("u_CloudBoxSize", Clouds::CloudRenderer::GetBoxSize());
		ColorShader.SetBool("u_CloudsEnabled", CloudsEnabled);
		ColorShader.SetBool("u_POM", POM);
		ColorShader.SetBool("u_HighQualityPOM", HighQualityPOM);
		ColorShader.SetBool("u_RTAO", RTAO);
		ColorShader.SetVector2f("u_Dimensions", glm::vec2(app.GetWidth(), app.GetHeight()));

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
		glBindTexture(GL_TEXTURE_2D, RoughReflections ? (DenoiseReflections ? ReflectionDenoised.GetTexture() : ReflectionTemporalFBO.GetTexture()) : ReflectionTraceFBO.GetTexture());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

		glActiveTexture(GL_TEXTURE12);
		glBindTexture(GL_TEXTURE_2D, CloudData);

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
		TemporalAAShader.SetInteger("u_PreviousPositionTexture", 3);
		TemporalAAShader.SetBool("u_Enabled", TAA);

		TemporalAAShader.SetMatrix4("u_PrevProjection", PreviousProjection);
		TemporalAAShader.SetMatrix4("u_PrevView", PreviousView);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetColorTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetPositionTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevTAAFBO.GetTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetPositionTexture());

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
			BloomRenderer::RenderBloom(BloomFBO, ColoredFBO.GetColorTexture(), ColoredFBO.GetPBRTexture());
		}

		// ---- Auto Exposure ----

		float ComputedExposure = 3.0f;

		if (AutoExposure)
		{
			#define SAMPLE_COUNT 9

			int iWidth = app.GetWidth();
			int iHeight = app.GetHeight();
			float AverageBrightness;

			int Pixels[SAMPLE_COUNT][2] =
			{
				static_cast<int>(iWidth * 0.50), static_cast<int>(iHeight * 0.50),
				static_cast<int>(iWidth * 0.25), static_cast<int>(iHeight * 0.50),
				static_cast<int>(iWidth * 0.75), static_cast<int>(iHeight * 0.50),
				static_cast<int>(iWidth * 0.50), static_cast<int>(iHeight * 0.25),
				static_cast<int>(iWidth * 0.50), static_cast<int>(iHeight * 0.75),
				static_cast<int>(iWidth * 0.25), static_cast<int>(iHeight * 0.25),
				static_cast<int>(iWidth * 0.25), static_cast<int>(iHeight * 0.75),
				static_cast<int>(iWidth * 0.75), static_cast<int>(iHeight * 0.25),
				static_cast<int>(iWidth * 0.75), static_cast<int>(iHeight * 0.75)
			};

			glm::vec3 kAveragedSamples;
			unsigned char Samples[SAMPLE_COUNT][4];

			for (unsigned int i = 0; i < SAMPLE_COUNT; i++)
			{
				glReadPixels(Pixels[i][0], Pixels[i][1], 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &Samples[i][0]);
				kAveragedSamples += (glm::vec3(Samples[i][0], Samples[i][1], Samples[i][2]) / 255.0f) / float(SAMPLE_COUNT);
			}

			AverageBrightness = glm::max(glm::max(kAveragedSamples.x, kAveragedSamples.y), kAveragedSamples.z);

			CameraExposure = 0.5 / AverageBrightness;
			CameraExposure = PrevCameraExposure + (CameraExposure - PrevCameraExposure) * 0.02;
			PrevCameraExposure = CameraExposure;

			ComputedExposure = CameraExposure * 3.0f;
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
		PostProcessingShader.SetInteger("u_RTAOTexture", 11);
		PostProcessingShader.SetInteger("u_NormalTexture", 12);
		PostProcessingShader.SetInteger("u_PBRTexture", 13);
		PostProcessingShader.SetInteger("u_CloudData", 15);
		PostProcessingShader.SetInteger("u_GodRaysStepCount", GodRaysStepCount);
		PostProcessingShader.SetVector3f("u_SunDirection", SunDirection);
		PostProcessingShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
		PostProcessingShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
		PostProcessingShader.SetVector2f("u_Dimensions", glm::vec2(PostProcessingFBO.GetWidth(), PostProcessingFBO.GetHeight()));
		PostProcessingShader.SetMatrix4("u_ProjectionMatrix", MainCamera.GetProjectionMatrix());
		PostProcessingShader.SetMatrix4("u_ViewMatrix", MainCamera.GetViewMatrix());
		PostProcessingShader.SetBool("u_SunIsStronger", StrongerLightDirection == SunDirection);
		PostProcessingShader.SetBool("u_LensFlare", LensFlare);
		PostProcessingShader.SetBool("u_GodRays", GodRays);
		PostProcessingShader.SetBool("u_SSAO", SSAO);
		PostProcessingShader.SetBool("u_Bloom", Bloom);
		PostProcessingShader.SetBool("u_SSGodRays", FakeGodRays);
		PostProcessingShader.SetBool("u_RTAO", RTAO);
		PostProcessingShader.SetBool("u_ExponentialFog", ExponentialFog);
		PostProcessingShader.SetBool("u_AutoExposure", AutoExposure);
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
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[0]);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[1]);

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[2]);

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[3]);
		//

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, ShadowFBO.GetTexture());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, VolumetricFBO.GetTexture());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D, RTAOTemporalFBO.GetTexture());

		glActiveTexture(GL_TEXTURE12);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetNormalTexture());

		glActiveTexture(GL_TEXTURE13);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetPBRTexture());

		glActiveTexture(GL_TEXTURE15);
		glBindTexture(GL_TEXTURE_2D, CloudData);

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

		// Particles

		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		if (RenderParticles)
		{
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			world->UpdateParticles(&MainCamera, InitialTraceFBO->GetPositionTexture(), ShadowFBO.GetTexture(), DiffuseDenoiseFBO.GetTexture(), SunDirection, MainCamera.GetPosition(), glm::vec2(app.GetWidth(), app.GetHeight()), DeltaTime);
			glDisable(GL_BLEND);
		}

		RendererUI.RenderQuad(glm::vec2(floor((float)app.GetWidth() / 2.0f), floor((float)app.GetHeight() / 2.0f)), &Crosshair, &OCamera);

		if (app.GetCurrentFrame() % 80 == 0)
		{
			world->GenerateDistanceField();
		}

		if (app.GetCurrentFrame() % 60 == 0)
		{
			world->m_ParticleEmitter.CleanUpList();
		}

		world->Update(&MainCamera);

		// Finish Frame
		glFinish();
		app.FinishFrame();
		
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glUseProgram(0);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		std::string title = "Voxel RT | "; title += BlockDatabase::GetBlockName(world->GetCurrentBlock()); title += "  ";
		GLClasses::DisplayFrameRate(app.GetWindow(), title);

		float CurrentTime = glfwGetTime();
		DeltaTime = CurrentTime - Frametime;
		Frametime = glfwGetTime();

		ModifiedWorld = false;
	}

	SaveWorld(world, world_name);
	delete world;

	return;
}
