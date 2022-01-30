// Main program render/gameplay pipeline code



// License

//////////////////////////////////////////////////////////////////////////////////////
//    																					
//    MIT License																		
//    																					
//    Copyright (c) 2021 Samuel Rasquinha												
//    																					
//    Permission is hereby granted, free of charge, to any person obtaining a copy		
//    of this software and associated documentation files (the "Software"), to deal		
//    in the Software without restriction, including without limitation the rights		
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell			
//    copies of the Software, and to permit persons to whom the Software is				
//    furnished to do so, subject to the following conditions:							
//    																					
//    The above copyright notice and this permission notice shall be included in all	
//    copies or substantial portions of the Software.									
//    																					
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR		
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,			
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE		
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER			
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,		
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE		
//    SOFTWARE.																			
//    																		            
//////////////////////////////////////////////////////////////////////////////////////


// License




// includes 

#include "Pipeline.h"
#include <chrono>
#include "ShaderManager.h"
#include "BlockDataSSBO.h"
#include "BlueNoiseDataSSBO.h"
#include "SoundManager.h"
#include "TAAJitter.h"
#include "VolumetricFloodFill.h"
#include "NBT/Importer.h"
#include "AnimatedTexture.h"

// includes


// static globals 
// these variables are used all over the file 
// yes, I know i'm a degenerate for using globals but whatever 

static VoxelRT::Player MainPlayer;

int VoxelRT_FloodFillDistanceLimit = 4; // Used in other files as an extern

static int DEBUGTraceLevel = 0;

static int DiffuseTraceLength = 42;
static int ReflectionTraceLength = 64;

static bool VSync = false;

static bool Fucktard = false;

static bool DiffuseDirectLightSampling = false;

static float DiffuseIndirectSuperSampleRes = 0.25f;

static bool JitterSceneForTAA = false;

static bool CloudReflections = true;

//static bool SmartUpscaleCloudTemporal = true;

static int LPVDebugState = 0;

static bool LightListDebug = false;

static bool HejlBurgessTonemap = false;

static float TextureDesatAmount = 0.1f;
static bool PreTemporalSpatialPass = true;
static float PurkingeEffectStrength = 0.0f;
static int SelectedColorGradingLUT = -1;
static bool ColorDither = true;
static float FilmGrainStrength = 0.0f;
static float ChromaticAberrationStrength = 0.0f;
static float ExposureMultiplier = 1.0f;

static float DiffractionResolution = 0.25f;
static float DiffractionStrength = 1.0f;
static int DiffractionKernel = 14;
static float DiffractionScaler = 1.0f;
static bool BlurDiffraction = false;

// pixel padding
static int PIXEL_PADDING = 20;

static bool TEMPORAL_SPEC = true;

static bool APPLY_PLAYER_SHADOW_FOR_GI = false;

static int RenderDistance = 400;

static bool REFLECT_PLAYER = false;

static bool DENOISE_REFLECTION_HIT_DATA = true;

static bool SVGF_LARGE_KERNEL = false;

static float ColorPhiBias = 3.325f;

//static bool ANTI_FLICKER = true;



static bool DenoiseSunShadows = true;


static bool AGGRESSIVE_DISOCCLUSION_HANDLING = true;

static bool CHECKERBOARD_SPP = true;
static bool CHECKERBOARD_SPEC_SPP = false;
static bool InferSpecularDetailSpatially = false;

static float GLOBAL_RESOLUTION_SCALE = 1.0f;

static bool ContrastAdaptiveSharpening = true;
static float CAS_SharpenAmount = 0.3 + 0.05f;

static bool UseEnvironmentBRDF = true;
static bool RemoveTiling = false;


static bool VXAO_CUTOFF = true;




// Clouds ->
static bool CloudsEnabled = true;
static float CloudCoverage = 1.1f;
static const bool CloudBayer = true;
static float CloudDetailScale = 1.1f;
static float CloudErosionWeightExponent = 0.56f;
static float CloudDetailFBMPower = 1.125f;
static bool CloudDetailWeightEnabled = false;
static bool ClampCloudTemporal = false;
static glm::vec2 CloudModifiers = glm::vec2(-0.950, 0.1250f); 
static bool CurlNoiseOffset = false;
static float CirrusScale = 1.8667f;
static float CirrusStrength = 0.11f;
static float CloudTimeScale = 1.0f;
static glm::ivec3 CloudStepCount = glm::ivec3(32, 8, 4);
static bool CloudCheckerStepCount = false;
static bool CloudLODLighting = true;
static bool CloudForceSupersample = true;
static float CloudForceSupersampleRes = 1.0f;
static float CloudResolution = 0.200f;
static bool CloudSpatialUpscale = true;
static bool CloudFinalCatmullromUpsample = false;
static float CloudAmbientDensityMultiplier = 1.2f;
static float CloudThiccness = 1250.0f;

// basic nebula settings ->
static float NebulaStrength = 1.0f;
static float NebulaGIStrength = 0.50f;
static bool NebulaCelestialColor = false;
static bool NebulaGIColor = true;

//

// Depth of field
static bool DOF = false;
static bool LargeKernelDOF = false;
static float DOFCOCScale = 2.5f;
static float DOFBlurScale = 1.0f;
static float DOFCAScale = 1.75f;
static float DOFResolution = 0.350f;
static float DOFTemporalDepthBlend = 0.875f;

static bool VXAO = true;
static bool WiderSVGF = false;
static bool DITHER_SPATIAL_UPSCALE = true;
static bool USE_NEW_SPECULAR_SPATIAL = true;
static bool USE_BLUE_NOISE_FOR_TRACING = true;

static bool PointVolumetricsToggled = false;
static float PointVolumetricStrength = 0.525f;
static bool ColoredPointVolumetrics = false;
static bool PointVolumetricsBayer = true;
static bool PointVolPerlinOD = false;
static bool DenoisePointVol = true;
static bool PointVolTriquadraticDensityInterpolation = false;
static bool PointVolGroundTruthColorInterpolation = false;

static float InitialTraceResolution = 1.0f;
static float GBufferResolution = 1.0f;
static float DiffuseTraceResolution = 0.250f; 

static float ShadowTraceResolution = 0.500f;
static float ReflectionTraceResolution = 0.250; 
static float ReflectionSuperSampleResolution = 0.250f;


static float SSAOResolution = 0.35f;
static float SSAOStrength = 0.325f;
static float RTAOResolution = 0.125f;
static float VolumetricResolution = 0.5f;

static float SunTick = 50.0f;
static float PreviousSunTick = 50.0f;
static float DiffuseLightIntensity = 1.2f;
static float LensFlareIntensity = 0.075f;
static float BloomQuality = 0.25f;
static bool BloomWide = true;
static float BloomStrength = 1.0f;
bool LensDirt = false;
float LensDirtStrength = 1.1f;
bool HQBloomUpscale = true;

static bool SoftShadows = true;
static float ShadowSupersampleRes = 0.5f;
static float ShadowFilterScale = 1.0f;


static bool SSSSS = false;
static float SSSSSStrength = 0.9f / 1.0f;


static int DiffuseSPP = 3; 
static int ReflectionSPP = 2;

// Alpha test : 
static bool ShouldAlphaTest = false;
static bool ShouldAlphaTestShadows = false;


static bool TAA = true;
const bool TAADepthWeight = true;
static float TAADepthWeightExp = 2.0f;
static bool FXAA = true;
static bool Bloom = true;
static bool DoDiffractionSpikes = false;
static bool USE_SVGF = true;
static bool DO_VARIANCE_SPATIAL = true;
static bool DO_SVGF_SPATIAL = true;
static bool DO_SVGF_TEMPORAL = true;


static bool BrutalFXAA = true;
//static bool DoSecondaryFXAA = true;


static bool GodRays = false;
static bool FakeGodRays = false;
static bool RoughReflections = true;
static bool DenoiseReflections = true;
static bool ReflectionFireflyRejection = true;
static bool AggressiveFireflyRejection = true;
static bool SmartReflectionClip = true;
static bool ReflectionNormalMapWeight = true;
static int ReflectionDenoisingRadiusBias = 0;
static bool ReflectionLPVGI = true;
static bool RefectionUseDecoupledGI = false;
static bool ReflectionHighQualityLPVGI = false;
static bool ReflectionRoughnessBias = true;
static bool ReflectionDenoiserDeviationHandling = true;
static float NormalMapWeightStrength = 0.750f;
static float ReflectionDenoiserScale = 1.0f;
static bool ReprojectReflectionsToScreenSpace = true;
static bool DeriveReflectionsFromDiffuseSH = false;
static bool TemporallyStabializeHitDistance = !false; // ;)
static bool ReflectionTemporalWeight = true;
static float RoughnessNormalWeightBiasStrength = 1.075f;
static bool AmplifyReflectionTransversalWeight = true;

static bool RenderParticles = true;

static bool LensFlare = false;
static bool SSAO = false;
static bool RTAO = false;
static bool POM = false;
static float POMHeight = 1.0f;
static bool DitherPOM = false;
static bool HighQualityPOM = false;

//static bool CheckerboardClouds = true;
static bool AmplifyNormalMap = true;


static bool DoSmallItemCube = true;
static bool SimpleLightingItemCube = false;
static int AntialiasItemCubeLevel = 4;
static int ItemCubeSpecularSampleBias = 0;


static int GodRaysStepCount = 12;
static float GodRaysStrength = 0.5f;

// Color modifiers
static float SunStrengthModifier = 0.850f;
static float MoonStrengthModifier = 1.0f;
static float GISunStrength = 1.0f;
static float GISkyStrength = 1.125f;

static bool AutoExposure = false;
static bool ExponentialFog = false;


static VoxelRT::World* world = nullptr;
static bool ModifiedWorld = false;
static VoxelRT::FPSCamera& MainCamera = MainPlayer.Camera;
static VoxelRT::OrthographicCamera OCamera(0.0f, 800.0f, 0.0f, 600.0f);

static float Frametime;
static float DeltaTime;

float VoxelRT_VolumeMultiplier = 1.0f;

static float DeltaSum = 0.0f;

static float PointVolumetricsScale = 0.2f; // 1/25 the pixels


static glm::vec3 CloudProjectionSunDir = glm::vec3(0.0f);



static glm::vec3 g_SunDirection;
static glm::vec3 g_MoonDirection;

static bool RandomDebugVar = false;


const bool DOWNSAMPLE_GBUFFERS = false;

static bool HighlightFocusedBlock = false;


static float CenterDepthSmooth = 1.0f;



std::array<glm::mat4, 6> GetMatrices(const glm::vec3& center) {

	std::array<glm::mat4, 6> view_matrices = {
		glm::lookAt(center, center + glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
		glm::lookAt(center, center + glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f,  0.0f,  1.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f,-1.0f, 0.0f), glm::vec3(0.0f,  0.0f, -1.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f, 0.0f,-1.0f), glm::vec3(0.0f, -1.0f,  0.0f))
	};

	return view_matrices;
}

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
			ImGui::Checkbox("Debug Variable", &RandomDebugVar);
			ImGui::Text("Player Position : %f, %f, %f", MainCamera.GetPosition().x, MainCamera.GetPosition().y, MainCamera.GetPosition().z);
			ImGui::Text("Camera Front : %f, %f, %f", MainCamera.GetFront().x, MainCamera.GetFront().y, MainCamera.GetFront().z);
			ImGui::Text("Player Grounded : %d", MainPlayer.m_isOnGround);
			ImGui::Text("Player Velocity : %f, %f, %f", MainPlayer.m_Velocity[0], MainPlayer.m_Velocity[1], MainPlayer.m_Velocity[2]);
			ImGui::Text("Sun Direction : %f, %f, %f", g_SunDirection[0], g_SunDirection[1], g_SunDirection[2]);
			ImGui::Text("Moon Direction : %f, %f, %f", g_MoonDirection[0], g_MoonDirection[1], g_MoonDirection[2]);
			ImGui::Text("Time : %f", glfwGetTime());
			ImGui::Text("Sun Tick : %f", SunTick);

			if (DOF)
				ImGui::Text("Temporal Center Depth (x1000) : %f", CenterDepthSmooth * 1000.0f);
			
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Checkbox("Highlight focused block?", &HighlightFocusedBlock);
			ImGui::SliderFloat("Mouse Sensitivity", &MainPlayer.Sensitivity, 0.025f, 1.0f);
			ImGui::SliderFloat("Player Speed", &MainPlayer.Speed, 0.025f, 0.3f); 

			ImGui::NewLine();
			if (ImGui::Button("Reset Speed/Sensitivity")) 
			{
				MainPlayer.Speed = 0.045f;
				MainPlayer.Sensitivity = 0.25f;
			}

			ImGui::NewLine();
			ImGui::Checkbox("VSync", &VSync);
			ImGui::NewLine();
			ImGui::NewLine();

			ImGui::SliderInt("TRACE DEBUG Level (0 : None, 1 : Indirect Diffuse, 2 : Specular, 3 : Direct", &DEBUGTraceLevel, 0, 3);

			ImGui::NewLine();
			ImGui::Checkbox("* Use accurate indirect light combining? (Option for stylization) (Accounts for metals)", &UseEnvironmentBRDF);
			ImGui::NewLine();

			// Sun/Moon ->
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Text("--- Sun/Moon Settings ---");
			ImGui::NewLine();
			ImGui::SliderFloat("Sun Time ", &SunTick, 0.1f, 256.0f);
			ImGui::NewLine();
			ImGui::SliderFloat("Sun Strength", &SunStrengthModifier, 0.01f, 1.0f);
			ImGui::SliderFloat("Moon Strength", &MoonStrengthModifier, 0.01f, 2.0f);
			ImGui::SliderFloat("GI Sun Strength", &GISunStrength, 0.01f, 2.5f);
			ImGui::SliderFloat("GI Sky Strength", &GISkyStrength, 0.01f, 2.5f);

			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Text("--- Item Cube Settings ---");
			ImGui::Checkbox("Render Item Cube", &DoSmallItemCube);
			ImGui::Checkbox("Simple Lighting", &SimpleLightingItemCube);
			ImGui::SliderInt("Specular Indirect Sample Bias (More = better reflections at the cost of performance)", &ItemCubeSpecularSampleBias, 0, 48);
			ImGui::SliderInt("AntiAliasing Item Cube Level", &AntialiasItemCubeLevel, 0, 16);


			// Initial hit
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Text("--- Initial hit Settings ---");
			ImGui::NewLine();
			ImGui::SliderFloat("Initial Trace Resolution", &InitialTraceResolution, 0.1f, 1.0f);
			ImGui::SliderFloat("GBuffer Generate Resolution", &GBufferResolution, 0.1f, 1.0f);
			ImGui::SliderInt("DF Trace Length.", &RenderDistance, 20, 500);
			ImGui::NewLine();


			// Indirect diffuse ->
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Text("--- Indirect Diffuse Settings ---");
			ImGui::NewLine();
			ImGui::Checkbox("WIP Light List Based Diffuse Direct Light Sampling", &DiffuseDirectLightSampling);
			ImGui::SliderFloat("Diffuse Trace Resolution ", &DiffuseTraceResolution, 0.1f, 1.25f);
			ImGui::SliderFloat("Diffuse Light Intensity ", &DiffuseLightIntensity, 0.05f, 1.25f);
			ImGui::NewLine();
			ImGui::SliderInt("Diffuse Trace LENGTH", &DiffuseTraceLength, 2, 96);
			ImGui::SliderInt("Diffuse Trace SPP", &DiffuseSPP, 1, 32);
			ImGui::NewLine();
			ImGui::Checkbox("Apply player shadow for global illumination?", &APPLY_PLAYER_SHADOW_FOR_GI);
			ImGui::Checkbox("CHECKERBOARD_DIFFUSE_SPP", &CHECKERBOARD_SPP);
			ImGui::SliderFloat("Diffuse Indirect Supersample Res", &DiffuseIndirectSuperSampleRes, 0.0f, 1.5f);
			ImGui::Checkbox("Pre Temporal Indirect Diffuse Spatial Pass?", &PreTemporalSpatialPass);
			ImGui::Checkbox("Use Blue noise for tracing?", &USE_BLUE_NOISE_FOR_TRACING);
			ImGui::NewLine();
			ImGui::Checkbox("Use SVGF? (Uses Atrous if disabled, SVGF recommended) ", &USE_SVGF);
			ImGui::SliderFloat("SVGF : Color Phi Bias", &ColorPhiBias, 0.5f, 6.0f);
			ImGui::Checkbox("AGGRESSIVE_DISOCCLUSION_HANDLING ", &AGGRESSIVE_DISOCCLUSION_HANDLING);
			ImGui::Checkbox("DO_SVGF_SPATIAL ", &DO_SVGF_SPATIAL);
			ImGui::Checkbox("SVGF_LARGE_KERNEL ", &SVGF_LARGE_KERNEL);
			ImGui::Checkbox("DO_SVGF_TEMPORAL ", &DO_SVGF_TEMPORAL);
			ImGui::Checkbox("DO_VARIANCE_SVGF_SPATIAL ", &DO_VARIANCE_SPATIAL);
			ImGui::Checkbox("WIDE_SVGF_SPATIAL ", &WiderSVGF);
			ImGui::NewLine();
			ImGui::NewLine();

			// Reflections ->
			ImGui::NewLine();
			ImGui::Text("--- Reflections/Indirect Specular Settings --");
			ImGui::NewLine();
			ImGui::SliderFloat("Reflection Trace Resolution ", &ReflectionTraceResolution, 0.1f, 1.25f);
			ImGui::NewLine();
			ImGui::SliderInt("Reflection Trace LENGTH", &ReflectionTraceLength, 2, 128);
			ImGui::SliderInt("Reflection Trace SPP", &ReflectionSPP, 1, 24);
			ImGui::NewLine();
			ImGui::SliderFloat("Reflection Super Sample Resolution", &ReflectionSuperSampleResolution, 0.05f, 1.5f);
			ImGui::Checkbox("Rough reflections?", &RoughReflections);
			
			ImGui::NewLine();
			ImGui::Checkbox("Temporally Filter Specular? (If turned off, increase SPP to stabialize)", &TEMPORAL_SPEC);
			ImGui::NewLine();
			ImGui::Checkbox("Use new reflection denoiser? (Preserves detail, might increase noise in some rare cases.)", &USE_NEW_SPECULAR_SPATIAL);
			ImGui::Checkbox("Denoise reflections?", &DenoiseReflections);
			ImGui::Checkbox("Strong Reflection Transversal Weight? (Improves contact hardening clarity, increases noise.)", &AmplifyReflectionTransversalWeight);
			ImGui::NewLine();
			ImGui::Checkbox("Denoise specular reprojection data? ", &DENOISE_REFLECTION_HIT_DATA);

			if (TEMPORAL_SPEC)
				ImGui::Checkbox("Temporally filter specular reprojection data? (gets rid of flickering artifacts from denoiser)", &TemporallyStabializeHitDistance);
			ImGui::SliderFloat("Reflection denoiser scale", &ReflectionDenoiserScale, 0.25f, 6.0f);
			ImGui::NewLine();

			if (TEMPORAL_SPEC)
				ImGui::Checkbox("Apply Reflection Denoiser Accumulation Weight?", &ReflectionTemporalWeight);
			
			ImGui::NewLine();
			ImGui::Checkbox("Derive reflections from diffuse spherical harmonic (derives when the material is too rough)", &DeriveReflectionsFromDiffuseSH);

			if (DenoiseReflections) {
				ImGui::Checkbox("Reflection Roughness Bias? (Increases reflection clarity by biasing the roughness value)", &ReflectionRoughnessBias);
				ImGui::Checkbox("Reflection Roughness-Based Distance Weight Clamping?", &ReflectionDenoiserDeviationHandling);
				
			}



			if (USE_NEW_SPECULAR_SPATIAL) {
				ImGui::SliderInt("Reflection Denoiser Radius Bias", &ReflectionDenoisingRadiusBias, -4, 4);
				ImGui::Checkbox("Normal map aware filtering? (Slightly more noise.)", &ReflectionNormalMapWeight);
				
				if (ReflectionNormalMapWeight)
					ImGui::SliderFloat("Normal map weight exponent", &NormalMapWeightStrength, 0.1f, 4.0f);
					ImGui::SliderFloat("Normal map - Roughness bias strength (Higher = lesser clarity with lesser noise, lower = more clarity with more noise", &RoughnessNormalWeightBiasStrength, 0.0f, 3.0f);
			}

			ImGui::NewLine();

			if (TEMPORAL_SPEC) {
				ImGui::Checkbox("Reflection radiance clipping? (Reduces ghosting, might cause more noise)", &SmartReflectionClip);
				ImGui::Checkbox("Filter fireflies in reflections?", &ReflectionFireflyRejection);

				if (ReflectionFireflyRejection) {
					ImGui::Checkbox("Aggressive Firefly Rejection?", &AggressiveFireflyRejection);
				}
			}
			//ImGui::Checkbox("Smart sharpen reflections?", &InferSpecularDetailSpatially);

			ImGui::NewLine();

			ImGui::Checkbox("Reflect player capsule?", &REFLECT_PLAYER);
			ImGui::Checkbox("Use screen space data for GI in reflections?", &ReprojectReflectionsToScreenSpace);
			ImGui::Checkbox("Approximate GI in reflections using LPV?", &ReflectionLPVGI);
			ImGui::Checkbox("Decouple GI (Decouples sky and point light gi)?", &RefectionUseDecoupledGI);

			if (ReflectionLPVGI) {
				ImGui::Checkbox("High Quality LPVGI in reflections?", &ReflectionHighQualityLPVGI);

			}

			ImGui::Checkbox("CHECKERBOARD_SPEC_SPP", &CHECKERBOARD_SPEC_SPP);
			ImGui::NewLine();

			// Shadows ->
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Text("--- Ray Traced Direct Shadows Settings --");
			ImGui::NewLine();
			ImGui::SliderFloat("Shadow Trace Resolution ", &ShadowTraceResolution, 0.1f, 1.25f);
			ImGui::SliderFloat("Shadow Supersample Resolution ", &ShadowSupersampleRes, 0.1f, 1.25f);

			if (DenoiseSunShadows)
				ImGui::SliderFloat("Shadow Filter Scale ", &ShadowFilterScale, 0.25f, 3.0f);

			ImGui::Checkbox("Contact Hardening Shadows?", &SoftShadows);
			ImGui::Checkbox("Denoise sun (or moon.) shadows?", &DenoiseSunShadows);
			ImGui::Checkbox("SCREEN SPACE subsurface scattering? (Or SSS)", &SSSSS);
			ImGui::SliderFloat("SCREEN SPACE SSS Strength", &SSSSSStrength, 0.0f, 6.0f);
			ImGui::NewLine();
			ImGui::NewLine();

			ImGui::Text("--- Misc Settings --");

			ImGui::SliderInt("Pixel Padding Amount", &PIXEL_PADDING, 0, 128);
			ImGui::SliderInt("God ray raymarch step count", &GodRaysStepCount, 8, 64);
			ImGui::Checkbox("BAYER 4x4 DITHER SPATIAL UPSCALE", &DITHER_SPATIAL_UPSCALE);
			ImGui::Checkbox("Particles?", &RenderParticles);
			ImGui::Checkbox("Amplify normal map?", &AmplifyNormalMap);
			ImGui::NewLine();
			ImGui::NewLine();

			ImGui::Text("--- Point Light Volumetrics Settings ---");
			ImGui::NewLine();
			ImGui::Checkbox("Point Volumetric Fog? ", &PointVolumetricsToggled);
			ImGui::Checkbox("Colored Fog?", &ColoredPointVolumetrics);
			ImGui::Checkbox("Use Bayer Matrix?", &PointVolumetricsBayer);
			ImGui::Checkbox("Use Perlin Noise FBM for Optical Depth?", &PointVolPerlinOD);
			ImGui::Checkbox("Use Ground Truth Color Interpolation? (EXPENSIVE! LEAVE OFF IF UNSURE.)", &PointVolGroundTruthColorInterpolation);
			ImGui::SliderFloat("Volumetric Render Resolution", &PointVolumetricsScale, 0.05f, 1.0f);
			ImGui::SliderFloat("Point Volumetrics Strength", &PointVolumetricStrength, 0.0f, 4.0f);
			ImGui::Checkbox("Denoise?", &DenoisePointVol);
			ImGui::Checkbox("Triquadratic Density Interpolation (Much slower, Gets rid of all bilinear interpolation artifacts)?", &PointVolTriquadraticDensityInterpolation);
			ImGui::SliderInt("Flood Fill Distance Limit", &VoxelRT_FloodFillDistanceLimit, 2, 8, "%d");
			ImGui::SliderInt("LPV Debug State (0 = OFF, 1 = LEVEL, 2 = LEVEL * COLOR, 3 = LEVEL AS INDIRECT, 4 = [LEVEL * COLOR] AS INDIRECT)", &LPVDebugState, 0, 4, "%d");
			
			if (ImGui::Button("Reset Propogation Settings (REPROPOGATES as well!)")) { 
				VoxelRT_FloodFillDistanceLimit = 4; 
				world->RepropogateLPV_();
			}



			//
			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(ImColor(128, 0, 0)));
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(ImColor(255, 0, 0)));

			if (ImGui::Button("Repropogate *ENTIRE* LPV")) {
				world->RepropogateLPV_();
			}

			ImGui::PopStyleColor(2);
			//


			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::NewLine();

			// Clouds ->

			ImGui::Text("--- Volumetric cloud Settings ---");
			ImGui::NewLine();
			ImGui::Checkbox("Volumetric Clouds?", &CloudsEnabled);
			ImGui::Checkbox("Volumetric Cloud Reflections? (Doesn't affect performance!)", &CloudReflections);
			ImGui::SliderFloat("Cloud Thiccness", &CloudThiccness, 100.0, 2000.0f);
			//ImGui::Checkbox("High Quality Clouds? (Doubles the ray march step count)", &CloudHighQuality);
			ImGui::Checkbox("Cloud Spatial Upscale", &CloudSpatialUpscale);
			ImGui::Checkbox("Curl Noise Offset?", &CurlNoiseOffset);
			ImGui::SliderInt("Raymarch Step Count", &CloudStepCount[0], 4, 64);
			ImGui::SliderInt("Lightmarch Step Count", &CloudStepCount[1], 2, 32);
			ImGui::SliderInt("Ambient Step Count", &CloudStepCount[2], 2, 32);
			ImGui::Checkbox("Cloud LOD Lighting? (Faster, uses LODs for acceleration)", &CloudLODLighting);

			//ImGui::Checkbox("Use smart checker cloud upscale (uses catmull rom if disabled)?", &SmartUpscaleCloudTemporal);
			ImGui::SliderFloat("Volumetric Cloud Density Multiplier", &CloudCoverage, 0.5f, 5.0f);
			ImGui::SliderFloat("Volumetric Cloud AMBIENT Density Multiplier", &CloudAmbientDensityMultiplier, 0.0f, 5.0f);
			ImGui::SliderFloat("Volumetric Cloud Resolution", &CloudResolution, 0.1f, 1.0f);
			ImGui::Checkbox("Cloud Force Supersample? (Force samples to given resolution)", &CloudForceSupersample);
			ImGui::SliderFloat("Cloud Force Supersample Resolution ", &CloudForceSupersampleRes, 0.1f, 1.0f);
			ImGui::SliderFloat2("Volumetric Cloud Modifiers", &CloudModifiers[0], -3.0f, 0.5);
			ImGui::SliderFloat("Volumetric Cloud Detail Scale", &CloudDetailScale, 0.0f, 2.0f);
			ImGui::SliderFloat("Detail Weight Exponent : ", &CloudErosionWeightExponent, 0.2f, 3.0f);
			ImGui::SliderFloat("Detail FBM Power : ", &CloudDetailFBMPower, 0.5f, 3.0f);

			//ImGui::Checkbox("Checkerboard clouds?", &CheckerboardClouds);
			ImGui::Checkbox("Clamp Cloud temporal?", &ClampCloudTemporal);
			ImGui::Checkbox("Cloud Checker Step Count?", &CloudCheckerStepCount);

			ImGui::SliderFloat("Cloud Time Scale", &CloudTimeScale, 0.2f, 3.0f);

			ImGui::SliderFloat("Fake Cirrus Scale", &CirrusScale, 0.5f, 3.0f);
			ImGui::SliderFloat("Fake Cirrus Strength", &CirrusStrength, 0.0f, 1.0f);

			ImGui::Checkbox("Final Cloud Catmullrom Upsample", &CloudFinalCatmullromUpsample);

			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::NewLine();


			// Post process ->


			ImGui::Text("--- Post Process Settings ---");
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Checkbox("Nebula affects moon color?", &NebulaCelestialColor);
			ImGui::Checkbox("Nebula affects GI colors?", &NebulaGIColor);

			ImGui::SliderFloat("Nebula Strength", &NebulaStrength, 0.0f, 6.0f);
			ImGui::SliderFloat("Nebula GI Strength", &NebulaGIStrength, 0.0f, 8.0f);
			ImGui::NewLine();

			if (USE_SVGF) {
				ImGui::Checkbox("Ray traced ambient occlusion? (RTAO)", &VXAO);
				ImGui::Checkbox("RTAO Cutoff? (Prevents artifacts if trace resolution is low, set it to false if you have a high diffuse trace resolution)", &VXAO_CUTOFF);
			}

			else {
				ImGui::Checkbox("Secondary RTAO (not recommended!)", &RTAO);
				ImGui::SliderFloat("Secondary RTAO Resolution", &RTAOResolution, 0.1f, 0.5f);
			}
			ImGui::NewLine();

			ImGui::NewLine();
			ImGui::Checkbox("Depth of field (DOF) ?", &DOF);
			ImGui::Checkbox("Large Kernel DOF?", &LargeKernelDOF);
			ImGui::SliderFloat("DOF COC Scale", &DOFCOCScale, 0.0f, 7.0f);
			ImGui::SliderFloat("DOF Blur Scale", &DOFBlurScale, 0.25f, 4.0f);
			ImGui::SliderFloat("DOF Chromatic Aberration Scale", &DOFCAScale, 0.0f, 4.0f);
			ImGui::SliderFloat("DOF Resolution", &DOFResolution, 0.1f, 1.0f);
			ImGui::NewLine();

			//ImGui::SliderFloat("DOF Temporal Depth Blend Weight", &DOFTemporalDepthBlend, 0.1f, 0.95f); ImGui::SameLine();
			//if (ImGui::Button("Reset")) { DOFTemporalDepthBlend = 0.875f;  }


			ImGui::NewLine();
			ImGui::Checkbox("Bloom?", &Bloom);
			ImGui::SliderFloat("Bloom Resolution ", &BloomQuality, 0.25f, 0.5f);
			ImGui::Checkbox("Wide Bloom?", &BloomWide);
			ImGui::SliderFloat("Bloom Strength", &BloomStrength, 0.0f, 1.25f);
			ImGui::Checkbox("HQ Bloom Upscale?", &HQBloomUpscale);
			ImGui::Checkbox("Lens Dirt?", &LensDirt);

			if (LensDirt) {
				ImGui::SliderFloat("Lens Dirt Strength", &LensDirtStrength, 0.001f, 4.0f);
			}
			ImGui::NewLine();
			ImGui::Checkbox("Bloom Diffraction Spikes? (WIP!) ", &DoDiffractionSpikes);
			ImGui::SliderFloat("Diffraction Resolution", &DiffractionResolution, 0.05f, 1.0f);
			ImGui::SliderFloat("Diffraction Strength", &DiffractionStrength, 0.25f, 10.0f);
			ImGui::SliderFloat("Diffraction Scale", &DiffractionScaler, 0.125f, 1.5f);
			ImGui::SliderInt("Diffraction Kernel Size", &DiffractionKernel, 4, 48);
			ImGui::Checkbox("Blur? ", &BlurDiffraction);
			ImGui::NewLine();
			ImGui::NewLine();

			ImGui::Checkbox("Auto Exposure (WIP!) ?", &AutoExposure);
			ImGui::SliderFloat("Exposure Multiplier", &ExposureMultiplier, 0.01f, 1.0f);
			ImGui::Checkbox("CAS (Contrast Adaptive Sharpening)", &ContrastAdaptiveSharpening);
			ImGui::SliderFloat("CAS SharpenAmount", &CAS_SharpenAmount, 0.0f, 0.99f);
			ImGui::SliderFloat("Desaturation Amount", &TextureDesatAmount, 0.0f, 1.0f);
			ImGui::Checkbox("Hejl Burgess Tonemap? (Uses ACES tonemap if disabled)", &HejlBurgessTonemap);
			ImGui::SliderFloat("Purkinje Effect Strength", &PurkingeEffectStrength, 0.0f, 1.0f);
			ImGui::NewLine();
			ImGui::SliderInt("Current Color Grading LUT (-1 = No grading)", &SelectedColorGradingLUT, -1, 9, "%d");
			ImGui::Checkbox("Emit Footstep Particles?", &MainPlayer.m_EmitFootstepParticles);
			ImGui::Checkbox("Color Dither", &ColorDither);
			ImGui::SliderFloat("Film Grain", &FilmGrainStrength, 0.0f, 0.04f);
			ImGui::SliderFloat("Chromatic Aberration (OFF if negative or zero)", &ChromaticAberrationStrength, -0.01f, 0.1f);
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Checkbox("Temporal Anti Aliasing", &TAA);
			//ImGui::Checkbox("TAA Depth Weight (Reduces ghosting)", &TAADepthWeight);

			if (TAADepthWeight)
				ImGui::SliderFloat("TAA Depth Weight Strength (Higher = Lesser ghosting, might cause higher aliasing)", &TAADepthWeightExp, 0.1f, 8.0f);

			ImGui::Checkbox("Jitter Projection Matrix For TAA? (small issues, right now :( ) ", &JitterSceneForTAA);
			ImGui::NewLine();
			ImGui::Checkbox("Fast Approximate Anti Aliasing", &FXAA);
			//ImGui::Checkbox("SECOND-PASS Fast Approximate Anti Aliasing", &DoSecondaryFXAA);
			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Checkbox("Lens Flare?", &LensFlare);
			ImGui::SliderFloat("Lens Flare Intensity ", &LensFlareIntensity, 0.05f, 1.5f);
			ImGui::Checkbox("(Implementation - 1) (WIP) God Rays? (Slower)", &GodRays);
			ImGui::Checkbox("(Implementation - 2) (WIP) God Rays? (faster, more crisp, Adjust the step count in the menu)", &FakeGodRays);
			ImGui::SliderFloat("God Ray Strength", &GodRaysStrength, 0.0f, 2.0f);
			ImGui::Checkbox("Exponential Fog?", &ExponentialFog);

			ImGui::NewLine();
			
			ImGui::NewLine();
			ImGui::NewLine();


			

			ImGui::NewLine();
			ImGui::NewLine();
			ImGui::Text("--- WIP and not-recommended stuff ---");
			ImGui::NewLine();

			ImGui::Checkbox("Screen Space Ambient Occlusion? (VXAO/RTAO recommended, ssao sucks.)", &SSAO);
			ImGui::SliderFloat("SSAO Render Resolution ", &SSAOResolution, 0.1f, 1.0f);
			ImGui::SliderFloat("SSAO Strength", &SSAOStrength, 0.1f, 2.0f);
			ImGui::Checkbox("Alpha Test? (WIP, has a few artifacts.) ", &ShouldAlphaTest);
			ImGui::Checkbox("Alpha Test Shadows? (WIP, has a few artifacts.)", &ShouldAlphaTestShadows);
			ImGui::Checkbox("POM? (VERY WORK IN PROGRESS, \
				The textures adapted from minecraft resource packs use a different parallax representation that needs to be handles)", &POM);
			ImGui::Checkbox("High Quality POM?", &HighQualityPOM);
			ImGui::Checkbox("Dither POM?", &DitherPOM);
			ImGui::SliderFloat("POM Depth", &POMHeight, 0.0f, 3.0f);
		} ImGui::End();

		if (LightListDebug) {


			if (ImGui::Begin("Light List Debug")) {

				glm::ivec3 block_loc = glm::ivec3(glm::floor(MainPlayer.m_Position));
				int cx = int(floor(float(block_loc.x) / float(16)));
				int cy = int(floor(float(block_loc.y) / float(16)));
				int cz = int(floor(float(block_loc.z) / float(16)));
				int OffsetArrayFetchLocation = (cz * 24 * 8) + (cy * 24) + cx;
				glm::ivec2 ChunkData = world->LightChunkOffsets[OffsetArrayFetchLocation];

				ImGui::Text("Player Original Position : %f  %f  %f", MainPlayer.m_Position.x, MainPlayer.m_Position.y, MainPlayer.m_Position.z);
				ImGui::Text("Player Floored Position : %d  %d  %d", block_loc.x, block_loc.y, block_loc.z);
				ImGui::Text("Player Chunk : %d  %d  %d", cx, cy, cz);
				ImGui::Text("Fetched Data : %d  %d", ChunkData.x, ChunkData.y);

				for (int i = 0; i < ChunkData.y; i++) {
					std::stringstream ss;
					glm::vec3 current_pos = glm::vec3(world->LightChunkData[ChunkData.x + i]);
					ImGui::Text("%f      %f      %f", current_pos[0], current_pos[1], current_pos[2]);
				}
			} ImGui::End();
		}

		if (ImGui::Begin("Other Settings and properties"))
		{
			static bool soundpack = true;
			std::string sndpackalt = soundpack ? "Alternate" : "Default";
			if (ImGui::Button(std::string("Switch Sound Pack To " + sndpackalt).c_str())) {
				soundpack = !soundpack;
				VoxelRT::SoundManager::SetPack(soundpack);
			}

			ImGui::NewLine();

			if (world) {

				std::string s = MainPlayer.m_isOnGround ? "Yes" : "No";
				std::string s1 = "Air, not grounded.";

				//if (MainPlayer.m_isOnGround) 
				{

					glm::ivec3 Idx = glm::ivec3(glm::floor(MainPlayer.m_Position));
					Idx.y -= 2;

					if (Idx.x > 0 && Idx.x < WORLD_SIZE_X - 1 && 
						Idx.y > 0 && Idx.y < WORLD_SIZE_Y - 1 && 
						Idx.z > 0 && Idx.z < WORLD_SIZE_Z - 1)
					{
						auto blockat = world->GetBlock((uint16_t)Idx.x, (uint16_t)Idx.y, (uint16_t)Idx.z);
						s1 = blockat.block > 0 ? VoxelRT::BlockDatabase::GetBlockName(blockat.block) : s1;
					}

				}

				ImGui::Text("Stood On Block : %s", s1.c_str());
				ImGui::SliderFloat("Volume", &VoxelRT_VolumeMultiplier, 0.0f, 3.5f);
			}

			ImGui::NewLine();
			ImGui::NewLine();

			ImGui::NewLine();
			ImGui::Text("--- Spec Presets ---\n");

			if (ImGui::Button("LOW") == true)
			{
				InitialTraceResolution = 1.0f;
				ShadowTraceResolution = 0.5f;
				DiffuseSPP = 3;
				ColorPhiBias = 2.0f;
				ReflectionTraceResolution = 0.25f;
				DiffuseTraceResolution = 0.25f;
				RTAO = false;
			}

			if (ImGui::Button("MEDIUM") == true)
			{
				InitialTraceResolution = 1.0f;
				ShadowTraceResolution = 0.75f;
				DiffuseSPP = 8;
				ColorPhiBias = 3.5f;
				ReflectionTraceResolution = 0.25f;
				DiffuseTraceResolution = 0.25f;
				RTAO = false;
			}

			if (ImGui::Button("HIGH") == true) {
				InitialTraceResolution = 1.0f;
				ShadowTraceResolution = 0.75f;
				DiffuseSPP = 2;
				ColorPhiBias = 3.25f;
				ReflectionTraceResolution = 0.5f;
				DiffuseTraceResolution = 0.5f;
				RTAO = false;
			}

			if (ImGui::Button("INSANE") == true) {
				InitialTraceResolution = 1.0f;
				ShadowTraceResolution = 1.0f;
				ReflectionTraceResolution = 0.5f;
				DiffuseTraceResolution = 0.5f;
				DiffuseSPP = 4;
				ColorPhiBias = 3.5f;
				RTAO = false;
			}

			ImGui::NewLine();
			ImGui::NewLine();

			

			if (ImGui::Button("Reset"))
			{
				MainPlayer.Sensitivity = 0.25f;
				MainPlayer.Speed = 0.045f;
				VoxelRT_VolumeMultiplier = 1.0f;
				VSync = false;
			}
		} ImGui::End();
	}

	void OnEvent(VoxelRT::Event e) override
	{
		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_SPACE)
		{
			MainPlayer.Jump();
		}

		if (e.type == VoxelRT::EventTypes::MouseMove && GetCursorLocked())
		{
			MainCamera.UpdateOnMouseMovement(GetCursorX(), GetCursorY());
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_LEFT && this->GetCursorLocked())
		{
			if (world)
			{
				ModifiedWorld = world->Raycast(0, MainCamera.GetPosition(), MainCamera.GetFront(), MainPlayer.m_Velocity, !MainPlayer.m_isOnGround, DeltaTime);
				
				if (ModifiedWorld)
					std::cout << "\n\nMODIFIED WORLD!\n\n";

			}
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_RIGHT && this->GetCursorLocked())
		{
			if (world)
			{
				ModifiedWorld = world->Raycast(1, MainCamera.GetPosition(), MainCamera.GetFront(), MainPlayer.m_Velocity, !MainPlayer.m_isOnGround, DeltaTime);
				
				if (ModifiedWorld)
					std::cout << "\n\nMODIFIED WORLD!\n\n";
			}
		}

		if (e.type == VoxelRT::EventTypes::MousePress && e.button == GLFW_MOUSE_BUTTON_MIDDLE && this->GetCursorLocked())
		{
			if (world)
			{
				world->Raycast(2, MainCamera.GetPosition(), MainCamera.GetFront(), MainPlayer.m_Velocity, !MainPlayer.m_isOnGround, DeltaTime);
			}
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_F1)
		{
			this->SetCursorLocked(!this->GetCursorLocked());
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_F10)
		{
			std::cout << "\n\n--REUPLOADED VOLUMETRIC VOLUME TO GPU--\n\n";
			VoxelRT::Volumetrics::Reupload();
			//world->RebufferLightList();
		}
		
		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_F11)
		{
			LightListDebug = !LightListDebug;
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

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_C)
		{
			MainPlayer.DisableCollisions = !MainPlayer.DisableCollisions;
		}

		if (e.type == VoxelRT::EventTypes::KeyPress && e.key == GLFW_KEY_ESCAPE)
		{
			VoxelRT::SaveWorld(world, world->m_Name);
			delete world;
			exit(0);
		}

		if (e.type == VoxelRT::EventTypes::WindowResize)
		{
			float pwx = this->GetWidth() + PIXEL_PADDING;
			float pwy = this->GetHeight() + PIXEL_PADDING;
			MainCamera.SetAspect((float)pwx / (float)pwy);
			OCamera.SetProjection(0.0f, pwx, 0.0f, pwy);
		}


		//world->RebufferLightList();
		
	}

};

GLClasses::Framebuffer InitialTraceFBO_1(16, 16, { {GL_R16F, GL_RED, GL_FLOAT, true, true}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false}, {GL_R32F, GL_RED, GL_FLOAT, true, true} }, false);
GLClasses::Framebuffer InitialTraceFBO_2(16, 16, { {GL_R16F, GL_RED, GL_FLOAT, true, true}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false}, {GL_R32F, GL_RED, GL_FLOAT, true, true} }, false);
GLClasses::Framebuffer HalfResGBuffer(16, 16, { {GL_R32F, GL_RED, GL_FLOAT, true, true}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false} }, false);
GLClasses::Framebuffer QuarterResGBuffer(16, 16, { {GL_R32F, GL_RED, GL_FLOAT, true, true}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false} }, false);

//GLClasses::Framebuffer GeneratedGBuffer(16, 16, { {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGB16F, GL_RED, GL_FLOAT, true, true}, {GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, false, false}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false} }, false);
GLClasses::Framebuffer GeneratedGBuffer(16, 16, { {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGB16F, GL_RED, GL_FLOAT, true, true}, {GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, false, false}, {GL_RED, GL_RED, GL_UNSIGNED_BYTE, false, false} }, false);

GLClasses::Framebuffer DiffuseRawTraceFBO(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT }, { GL_RG, GL_RG, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer DiffuseTemporalFBO1(16, 16, {{ GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT }, { GL_RGB16F, GL_RGB, GL_FLOAT } , { GL_RG, GL_RG, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer DiffuseTemporalFBO2(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT }, { GL_RGB16F, GL_RGB, GL_FLOAT } , { GL_RG, GL_RG, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer DiffusePreTemporal_SpatialFBO(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT } , { GL_R16F, GL_RED, GL_FLOAT } , { GL_RG, GL_RG, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer DiffuseDenoiseFBO(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT } , { GL_R16F, GL_RED, GL_FLOAT } , { GL_RG, GL_RG, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer DiffuseDenoisedFBO2(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT } , { GL_R16F, GL_RED, GL_FLOAT } , { GL_RG, GL_RG, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer VarianceFBO(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_RG16F, GL_RG, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT } }, false);



GLClasses::Framebuffer VolumetricsCompute(16, 16, { { GL_RGB16F, GL_RGB, GL_FLOAT } }, false, false);
GLClasses::Framebuffer VolumetricsComputeBlurred(16, 16, { { GL_RGB16F, GL_RGB, GL_FLOAT } }, false, false);

//GLClasses::Framebuffer PostProcessingFBO_History[4] = {
//	GLClasses::Framebuffer(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false),
//	GLClasses::Framebuffer(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false),
//	GLClasses::Framebuffer(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false),
//	GLClasses::Framebuffer(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false) };

GLClasses::Framebuffer PostProcessingFBO(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false);
GLClasses::Framebuffer DOFFBO(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false);
GLClasses::Framebuffer DOFFBOInput(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false);
GLClasses::Framebuffer TonemappedFBO(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT }, false);
GLClasses::Framebuffer FXAA_FBO(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false);

//GLClasses::Framebuffer AntiFlickerFBO(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, false);



// lut
// (x % y) + (t % y)  
const glm::ivec3 HistoryLUT[4] = {
	glm::ivec3(1, 2, 3),
	glm::ivec3(2, 3, 0),
	glm::ivec3(3, 0, 1),
	glm::ivec3(0, 1, 2)
};


GLClasses::Framebuffer ReflectionTraceFBO_1(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT }, { GL_RED, GL_RED, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer ReflectionTraceFBO_2(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT }, { GL_RED, GL_RED, GL_UNSIGNED_BYTE } }, false);
GLClasses::Framebuffer ReflectionTemporalFBO_1(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT } }, false);
GLClasses::Framebuffer ReflectionTemporalFBO_2(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT }, { GL_R16F, GL_RED, GL_FLOAT } }, false);
GLClasses::Framebuffer ReflectionDenoised_1(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT } }, false);
GLClasses::Framebuffer ReflectionDenoised_2(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT } }, false);
GLClasses::Framebuffer ReflectionHitDataDenoised(16, 16, { { GL_R16F, GL_RED, GL_FLOAT } }, false);



// shadow buffers
GLClasses::Framebuffer ShadowRawTrace(16, 16, { { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, { GL_R16F, GL_RED, GL_FLOAT } }, false),
ShadowTemporalFBO_1(16, 16, { { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, { GL_R16F, GL_RED, GL_FLOAT } }, false), ShadowTemporalFBO_2(16, 16,{ { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, { GL_R16F, GL_RED, GL_FLOAT } }, false),
ShadowFiltered(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, false), ShadowSSS(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, false), ShadowSSS2(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, false);



GLClasses::Framebuffer BloomCombined(16, 16, { { GL_RGB16F, GL_RGB, GL_FLOAT } }, false);
GLClasses::Framebuffer DiffractionSpikes(16, 16, { { GL_RGB16F, GL_RGB, GL_FLOAT } }, false);
GLClasses::Framebuffer DiffractionSpikesDenoised(16, 16, { { GL_RGB16F, GL_RGB, GL_FLOAT } }, false);



void VoxelRT::MainPipeline::StartPipeline()
{
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	std::vector<glm::ivec3> LightLocations;

	RayTracerApp app;
	app.Initialize();
	std::cout << "\n\n--------CONTROLS-------- \n\tF1 -> Lock/Unlock mouse\
		\n\tW, S, A, D->Move\
		\n\tSpace->Jump / Fly\
		\n\tShift->Accelerate down\
		\n\tQ / E->Change current block\
		\n\tF->Toggle freefly\
		\n\tTAB->Boost Player Speed\
		\n\tC->Toggle Collisions(only toggles IF on freefly mode)\
		\n\tESC->Save and quit\
		\n\tV->Toggle VSync\
		\n\tF2->Recompile shaders\
		\n\tImGui Windows ->\
		\n\tWindow 1 : Various Other Settings(Resolution settings)\
		\n\tWindow 2 : Player sensitivity, speed and sound options\n\n ";
	
	VoxelRT::BlockDatabase::Initialize();

	bool gen_type = 0;
	int create_type = 0;

	std::string world_name;

	world = new VoxelRT::World();

	do {
		std::cout << "\nEnter the name of your world : ";
		std::cin >> world_name;
	} while (!VoxelRT::FilenameValid(world_name));

	world->m_Name = world_name;

	if (!LoadWorld(world, world_name, LightLocations))
	{
		std::cout << "\nWhat would you like to create your world with? (0 : TERRAIN GENERATOR, 1 : IMPORT MINECRAFT WORLD) : ";
		std::cin >> create_type;
		std::cout << "\n\n";

		if (create_type == 0) {

			std::cout << "\nWhat type of world would you like to generate? (FLAT = 0, PLAINS = 1) : ";
			std::cin >> gen_type;
			std::cout << "\n\n";

			GenerateWorld(world, gen_type);
		}

		else if (create_type == 1) {
			std::string MinecraftWorldPath = "\\";

			do
			{
				std::cout << "\n\nEnter the path of a valid readable/accessible minecraft *SAVE REGION* directory : ";
				std::cin.ignore();
				std::getline(std::cin, MinecraftWorldPath);
			} while (!VoxelRT::DirectoryValid(MinecraftWorldPath));

			glm::vec3 ImportOrigin = glm::vec3(0.0f);
			std::cout << "\nEnter the import origin (X Y Z) : ";
			std::cin >> ImportOrigin.x;
			std::cin >> ImportOrigin.y;
			std::cin >> ImportOrigin.z;
			VoxelRT::MCWorldImporter::ImportWorld(MinecraftWorldPath, &world->m_WorldData, ImportOrigin);
		}
	}


	int HardwareProfile = 0;

	std::cout << "\nHardware Spec? (0 -> Low, 1 -> Medium, 2 -> High, 3 -> Insane) : ";
	std::cin >> HardwareProfile;

	std::cout << "\n\n\n";

	if (HardwareProfile == 0)
	{
		// Defaults.
	}

	if (HardwareProfile == 1)
	{
		InitialTraceResolution = 1.0f;
		ShadowTraceResolution = 0.75f;
		DiffuseSPP = 6;
		ColorPhiBias = 3.5f;
		ReflectionTraceResolution = 0.25f;
		DiffuseTraceResolution = 0.25f;
		RTAO = false;
	}

	if (HardwareProfile == 2)
	{
		InitialTraceResolution = 1.0f;
		ShadowTraceResolution = 0.75f;
		DiffuseSPP = 2;
		ColorPhiBias = 3.25f;
		SVGF_LARGE_KERNEL = true;
		WiderSVGF = true;
		ReflectionTraceResolution = 0.350f;
		DiffuseTraceResolution = 0.5f;
		RTAO = false;
	}

	if (HardwareProfile == 3)
	{
		InitialTraceResolution = 1.0f;
		ShadowTraceResolution = 1.0f;
		ReflectionTraceResolution = 0.5f;
		DiffuseTraceResolution = 0.5f;
		DiffuseSPP = 4;
		ColorPhiBias = 3.5f;
		SVGF_LARGE_KERNEL = true;
		WiderSVGF = true;
		RTAO = false;
	}

	// Initialize world, df generator etc 
	world->Buffer();
	world->InitializeDistanceGenerator();
	world->GenerateDistanceField();

	// Initialize sound engine

	std::cout << "\n\n";
	std::cout << "Initializing Sound Engine..\n";
	SoundManager::InitializeSoundManager();
	std::cout << "\nInitialized Sound Engine.";
	std::cout << "\n\nLoading & Playing Sounds..\n";
	SoundManager::LoadSounds();
	std::cout << "\nLoaded Sounds..";
	std::cout << "\n\n";


	// Create and compile shaders 
	ShaderManager::CreateShaders();


	VoxelRT::Renderer2D RendererUI;
	GLClasses::VertexBuffer VBO;
	GLClasses::VertexArray VAO;
	GLClasses::Shader& InitialTraceShader = ShaderManager::GetShader("INITIAL_TRACE");
	GLClasses::Shader& FinalShader = ShaderManager::GetShader("FINAL_SHADER");
	GLClasses::Shader& DiffuseTraceShader = ShaderManager::GetShader("DIFFUSE_TRACE");
	GLClasses::Shader& MainTemporalFilter = ShaderManager::GetShader("MAIN_TEMPORAL_FILER");
	GLClasses::Shader& ShadowTemporalFilter = ShaderManager::GetShader("SHADOW_TEMPORAL");
	GLClasses::Shader& ColorShader = ShaderManager::GetShader("COLOR_SHADER");
	GLClasses::Shader& PostProcessingShader = ShaderManager::GetShader("POST_PROCESS");
	GLClasses::Shader& TemporalAAShader = ShaderManager::GetShader("TEMPORAL_AA");
	GLClasses::Shader& ShadowTraceShader = ShaderManager::GetShader("SHADOW_TRACE");
	GLClasses::Shader& ReflectionTraceShader = ShaderManager::GetShader("REFLECTION_TRACE");
	GLClasses::Shader& SSAOShader = ShaderManager::GetShader("SSAO");
	GLClasses::Shader& SSAO_Blur = ShaderManager::GetShader("SSAO_BLUR");
	GLClasses::Shader& SimpleDownsample = ShaderManager::GetShader("SIMPLE_DOWNSAMPLE");
	GLClasses::Shader& LumaAverager = ShaderManager::GetShader("LUMA_AVERAGER");
	GLClasses::Shader& VolumetricScattering = ShaderManager::GetShader("VOLUMETRIC_SCATTERING");
	GLClasses::Shader& BilateralBlur = ShaderManager::GetShader("BILATERAL_BLUR");
	GLClasses::Shader& ReflectionDenoiserOld = ShaderManager::GetShader("REFLECTION_DENOISER");
	GLClasses::Shader& ReflectionDenoiserNew = ShaderManager::GetShader("REFLECTION_DENOISER_NEW");
	GLClasses::Shader& RTAOShader = ShaderManager::GetShader("RTAO");
	GLClasses::Shader& GaussianSpatialFilter = ShaderManager::GetShader("GAUSSIAN_SPATIAL_FILTER");
	GLClasses::Shader& AtrousSpatialFilter = ShaderManager::GetShader("ATROUS_SPATIAL_FILTER");
	GLClasses::Shader& SpatialInitial = ShaderManager::GetShader("SPATIAL_INITIAL");
	GLClasses::Shader& SpecularTemporalFilter = ShaderManager::GetShader("SPECULAR_TEMPORAL");
	GLClasses::Shader& DefaultCheckerboardReconstructor = ShaderManager::GetShader("CHECKER_RECONSTRUCT");
	GLClasses::Shader& SpecularCheckerboardReconstructor = ShaderManager::GetShader("SPECULAR_CHECKER_RECONSTRUCT");
	GLClasses::Shader& ShadowFilter = ShaderManager::GetShader("SHADOW_FILTER");
	GLClasses::Shader& VarianceEstimator = ShaderManager::GetShader("VARIANCE_ESTIMATOR");
	GLClasses::Shader& PointVolumetrics = ShaderManager::GetShader("VOLUMETRICS_COMPUTE");
	GLClasses::Shader& VolumetricsDenoiser = ShaderManager::GetShader("VOLUMETRICS_DENOISER");
	GLClasses::Shader& Gaussian9TapOptimized = ShaderManager::GetShader("GAUSSIAN_9TAP_OPTIMIZED");
	GLClasses::Shader& Gaussian5TapOptimized = ShaderManager::GetShader("GAUSSIAN_5TAP_OPTIMIZED");
	GLClasses::Shader& BilateralHitDist_1 = ShaderManager::GetShader("BILATERAL_HITDIST1");
	GLClasses::Shader& BilateralHitDist_2 = ShaderManager::GetShader("BILATERAL_HITDIST2");
	GLClasses::Shader& SSSBlur = ShaderManager::GetShader("SSSBLUR");
	GLClasses::Shader& SSSBlurPoisson = ShaderManager::GetShader("SSSBLUR_POISSON");
	GLClasses::Shader& GBuffer4xDownsample = ShaderManager::GetShader("GBUFFER_DOWNSAMPLE");
	
	// wip.
	GLClasses::Shader& SVGF_Temporal = ShaderManager::GetShader("SVGF_TEMPORAL");
	GLClasses::Shader& SVGF_Spatial = ShaderManager::GetShader("SVGF_SPATIAL");
	GLClasses::Shader& SVGF_Variance = ShaderManager::GetShader("SVGF_VARIANCE");


	GLClasses::Shader& CAS_Shader = ShaderManager::GetShader("CONTRAST_ADAPTIVE_SHARPENING");
	GLClasses::Shader& AntiFlickerShader = ShaderManager::GetShader("ANTI_FLICKER");
	GLClasses::Shader& CubeItemRenderer = ShaderManager::GetShader("CUBE_ITEM_RENDERER");
	GLClasses::Shader& FXAA_Secondary = ShaderManager::GetShader("FXAA_SECONDARY");
	GLClasses::Shader& DOFShader = ShaderManager::GetShader("DOF");
	GLClasses::Shader& BicubicDownsampler = ShaderManager::GetShader("BICUBIC_DOWNSAMPLE");
	GLClasses::Shader& Tonemapper = ShaderManager::GetShader("TONEMAPPER");


	GLClasses::Shader& GenerateGBuffer = ShaderManager::GetShader("GBUFFER_GENERATE");
	GLClasses::Shader& CombineBloom = ShaderManager::GetShader("BLOOM_COMBINE");
	GLClasses::Shader& DiffractionSpikesShader = ShaderManager::GetShader("DIFFRACTION_SPIKES");

	GLClasses::TextureArray BlueNoise;


	
	VoxelRT::ColorPassFBO ColoredFBO;

	GLClasses::Framebuffer FXAA_Final(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, true);

	GLClasses::Framebuffer TAAFBO1(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, true);
	GLClasses::Framebuffer TAAFBO2(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT }, true);
	GLClasses::Framebuffer DownsampledFBO(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT }, false);
	GLClasses::Framebuffer AverageLumaFBO(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT }, false);
	GLClasses::Framebuffer VolumetricFBO(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true);
	GLClasses::Framebuffer BlurredVolumetricFBO(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true);
	GLClasses::Framebuffer RTAO_FBO(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true), RTAO_TemporalFBO_1(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true), RTAO_TemporalFBO_2(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true);
	VoxelRT::BloomFBO BloomFBO(16, 16);
	VoxelRT::BloomFBO BloomFBOAlternate(16, 16);
	GLClasses::Framebuffer SSAOFBO(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true);
	GLClasses::Framebuffer SSAOBlurred(16, 16, { GL_RED, GL_RED, GL_UNSIGNED_BYTE }, true);

	glm::mat4 CurrentProjection, CurrentView;
	glm::mat4 PreviousProjection, PreviousView;
	glm::mat4 ShadowProjection, ShadowView;
	glm::mat4 ReflectionProjection, ReflectionView;
	glm::vec3 CurrentPosition, PreviousPosition;

	VoxelRT::AtmosphereRenderMap SkymapMain(64); 
	VoxelRT::AtmosphereRenderMap SkymapSecondary(24); // 8 * 3
	VoxelRT::AtmosphereRenderMap SkymapSecondary_2(16); 
	VoxelRT::AtmosphereRenderer AtmosphereRenderer;

	BlueNoiseDataSSBO BlueNoise_SSBO;

	GLClasses::Texture Crosshair;
	GLClasses::Texture BluenoiseTexture;
	GLClasses::Texture BluenoiseHighResTexture;
	GLClasses::Texture PlayerSprite;
	GLClasses::Texture ColorGradingLUT;
	GLClasses::Texture BlueNoiseLowResolution;
	GLClasses::Texture LensDirtTexture;
	GLClasses::CubeTextureMap NightSkyMap;


	NightSkyMap.CreateCubeTextureMap({
			"Res/nebulae/right.png",
			"Res/nebulae/left.png",
			"Res/nebulae/top.png",
			"Res/nebulae/bottom.png",
			"Res/nebulae/front.png",
			"Res/nebulae/back.png"
	}, true);

	GLClasses::CubeTextureMap NightSkyMapLowRes;


	NightSkyMapLowRes.CreateCubeTextureMap({
			"Res/nebulae_lowres/right.png",
			"Res/nebulae_lowres/left.png",
			"Res/nebulae_lowres/top.png",
			"Res/nebulae_lowres/bottom.png",
			"Res/nebulae_lowres/front.png",
			"Res/nebulae_lowres/back.png"
		}, true);


	GLClasses::CubeTextureMap RandomHDRI;
	RandomHDRI.CreateCubeTextureMap({
			"Res/AssHdr/px.png",
			"Res/AssHdr/nx.png",
			"Res/AssHdr/py.png",
			"Res/AssHdr/ny.png",
			"Res/AssHdr/pz.png",
			"Res/AssHdr/nz.png"
		}, true);


	GLClasses::CubeTextureMap RandomHDRIDiffuse;
	RandomHDRIDiffuse.CreateCubeTextureMap({
			"Res/AssHdrDiffuse/px.png",
			"Res/AssHdrDiffuse/nx.png",
			"Res/AssHdrDiffuse/py.png",
			"Res/AssHdrDiffuse/ny.png",
			"Res/AssHdrDiffuse/pz.png",
			"Res/AssHdrDiffuse/nz.png"
		}, true);



	Crosshair.CreateTexture("Res/Misc/crosshair.png", false, false);
	ColorGradingLUT.CreateTexture("Res/Misc/colorluts.png", false);
	BluenoiseTexture.CreateTexture("Res/Misc/blue_noise.png", false);
	BluenoiseHighResTexture.CreateTexture("Res/Misc/BluenoiseHighRes.png", false);
	PlayerSprite.CreateTexture("Res/Misc/player.png", false, true);
	//BlueNoise64x.CreateTexture("Res/Misc/blue_noise64x.png", false, false);
	BlueNoiseLowResolution.CreateTexture("Res/Misc/blue_noise32x.png", false, false);
	LensDirtTexture.CreateTexture("Res/Misc/lensdirt.png", true, true);

	BlockDataSSBO BlockDataStorageBuffer;
	BlockDataStorageBuffer.CreateBuffers();
	
	BlueNoise.CreateArray({
		"Res/Misc/BL_0.png",
		"Res/Misc/BL_1.png",
		"Res/Misc/BL_2.png",
		"Res/Misc/BL_3.png"
		}, { 256, 256 }, true, false);


	AnimatedTexture LavaAlbedo;
	LavaAlbedo.Create("Res/Block/Lava/Frames/Albedo", 256, 7);
	AnimatedTexture LavaNormals;
	LavaNormals.Create("Res/Block/Lava/Frames/Normal", 256, 7);


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

	MainCamera.SetPosition(glm::vec3(WORLD_SIZE_X / 2, 75, WORLD_SIZE_Z / 2));

	BloomRenderer::Initialize();
	AverageLumaFBO.SetSize(1, 1);

	Clouds::CloudRenderer::Initialize();

	auto* InitialTraceFBO = &InitialTraceFBO_1;
	auto* InitialTraceFBOPrev = &InitialTraceFBO_2;

	glm::vec3 StrongerLightDirection;

	GLfloat PreviousLuma = 3.0f;

	float CameraExposure = 1.0f;
	float PrevCameraExposure = 1.0f;


	Frametime = glfwGetTime();

	bool UpdatePlayerCollision = true;

	GenerateJitterStuff();


	// Volumetricssss 

	world->InitializeLightList();

	Volumetrics::CreateVolume(world, BlockDataStorageBuffer.GetSSBO(), BlockDatabase::GetTextureArray());
	for (auto& e : LightLocations) {
		uint8_t block_at = world->GetBlock(e).block;
		Volumetrics::AddLightToVolume(e, block_at);
		world->InsertToLightList(e);
	}

	for (int i = 0; i < 3; i++) {
		Volumetrics::PropogateVolume();
	}


	/////////////

	// Auto exposure!

	const float ZeroFloat = 0.0f;
	GLuint AutoExposureSSBO = 0, AutoExposurePDF;
	GLuint Zero256UINT[256];
	memset(&Zero256UINT, 0, 256 * sizeof(GLuint));

	glGenBuffers(1, &AutoExposureSSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, AutoExposureSSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float), &ZeroFloat, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	glGenBuffers(1, &AutoExposurePDF);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, AutoExposurePDF);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLuint) * 256, &Zero256UINT, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);


	// DOF 

	GLuint DOFSSBO = 0;
	float DOFDATA = 0.;
	glGenBuffers(1, &DOFSSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, DOFSSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * 1, &DOFDATA, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);




	
	// Create Compute Shaders!
	GLClasses::ComputeShader AutoExposureComputePDF;
	AutoExposureComputePDF.CreateComputeShader("Core/Shaders/CalculateAverageLuminance.comp");
	AutoExposureComputePDF.Compile();

	GLClasses::ComputeShader ComputeAutoExposure;
	ComputeAutoExposure.CreateComputeShader("Core/Shaders/ComputeExposure.comp");
	ComputeAutoExposure.Compile();

	GLClasses::ComputeShader WriteCenterDepth;
	WriteCenterDepth.CreateComputeShader("Core/Shaders/WriteCenterDepth.comp");
	WriteCenterDepth.Compile();

	////////////

	// Equi rectangular projected cloud map ->



	GLClasses::Framebuffer CloudProjection(512 * 3, 512 * 3, { GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE }, true);
	CloudProjection.CreateFramebuffer();

	CloudProjection.Bind();
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	CloudProjection.Unbind();

	// 

	glm::ivec4 FocusedOnBlock = glm::ivec4(-1);

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		ShadowSupersampleRes = glm::max(ShadowSupersampleRes, ShadowTraceResolution);

		glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);

		auto& ReflectionDenoiser = USE_NEW_SPECULAR_SPATIAL ? ReflectionDenoiserNew : ReflectionDenoiserOld;

		// Player update flag
		if (glfwGetWindowAttrib(app.GetWindow(), GLFW_FOCUSED) == 0) {
			UpdatePlayerCollision = false;
		}

		else {
			UpdatePlayerCollision = true;
		}

		if (app.GetCurrentFrame() < 5) {
			UpdatePlayerCollision = true;
		}


		// Sound update ->

		SoundManager::UpdatePosition(MainCamera.GetFront(), MainCamera.GetPosition(), MainCamera.GetUp());



		// view bob


		// Jitter



		// Tick the sun and moon
		float time_angle = SunTick * 2.0f;
		glm::mat4 sun_rotation_matrix;

		glm::vec3 SunDirection;
		glm::vec3 MoonDirection;

		sun_rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(time_angle), glm::vec3(0.0f, 0.0f, 1.0f));
		SunDirection = glm::vec3(sun_rotation_matrix * glm::vec4(1.0f));
		MoonDirection = glm::vec3(-SunDirection.x, -SunDirection.y, SunDirection.z);
		StrongerLightDirection = -SunDirection.y < 0.01f ? SunDirection : MoonDirection;

		SunDirection = glm::normalize(SunDirection);
		MoonDirection = glm::normalize(MoonDirection);
		g_SunDirection = SunDirection;
		g_MoonDirection = MoonDirection;
		StrongerLightDirection = glm::normalize(StrongerLightDirection);


		glfwSwapInterval((int)VSync);
		
		float PADDED_WIDTH = app.GetWidth() + PIXEL_PADDING;
		float PADDED_HEIGHT = app.GetHeight() + PIXEL_PADDING;
		float TRUE_PADDED_WIDTH = app.GetWidth() + PIXEL_PADDING;
		float TRUE_PADDED_HEIGHT = app.GetHeight() + PIXEL_PADDING;
		PADDED_WIDTH *= GLOBAL_RESOLUTION_SCALE;
		PADDED_HEIGHT *= GLOBAL_RESOLUTION_SCALE;

		
		// Resize the framebuffers
		{
			
			// Without padding!
			FXAA_Final.SetSize(app.GetWidth(), app.GetHeight());



			InitialTraceFBO_1.SetSize(floor(TRUE_PADDED_WIDTH * InitialTraceResolution), floor(TRUE_PADDED_HEIGHT * InitialTraceResolution));
			InitialTraceFBO_2.SetSize(floor(TRUE_PADDED_WIDTH * InitialTraceResolution), floor(TRUE_PADDED_HEIGHT * InitialTraceResolution));
			GeneratedGBuffer.SetSize(floor(TRUE_PADDED_WIDTH * GBufferResolution), floor(TRUE_PADDED_HEIGHT * GBufferResolution));
			
			if (DOWNSAMPLE_GBUFFERS) {
				HalfResGBuffer.SetSize(floor(TRUE_PADDED_WIDTH * InitialTraceResolution * 0.5f), floor(TRUE_PADDED_HEIGHT * InitialTraceResolution * 0.5f));
				QuarterResGBuffer.SetSize(floor(TRUE_PADDED_WIDTH * InitialTraceResolution * 0.25f), floor(TRUE_PADDED_HEIGHT * InitialTraceResolution * 0.25f));
			}

			VolumetricsCompute.SetSize(floor(PADDED_WIDTH * PointVolumetricsScale), floor(PADDED_HEIGHT * PointVolumetricsScale));
			VolumetricsComputeBlurred.SetSize(floor(PADDED_WIDTH * PointVolumetricsScale), floor(PADDED_HEIGHT * PointVolumetricsScale));
			
			DiffuseRawTraceFBO.SetSize(PADDED_WIDTH * DiffuseTraceResolution, PADDED_HEIGHT * DiffuseTraceResolution);
			DiffusePreTemporal_SpatialFBO.SetSize(PADDED_WIDTH* DiffuseTraceResolution, PADDED_HEIGHT* DiffuseTraceResolution);
			
			DiffuseIndirectSuperSampleRes = glm::max(DiffuseIndirectSuperSampleRes, DiffuseTraceResolution);
			DiffuseTemporalFBO1.SetSize(PADDED_WIDTH * DiffuseIndirectSuperSampleRes, PADDED_HEIGHT * DiffuseIndirectSuperSampleRes);
			DiffuseTemporalFBO2.SetSize(PADDED_WIDTH * DiffuseIndirectSuperSampleRes, PADDED_HEIGHT * DiffuseIndirectSuperSampleRes);
			DiffuseDenoiseFBO.SetSize(PADDED_WIDTH * DiffuseIndirectSuperSampleRes, PADDED_HEIGHT * DiffuseIndirectSuperSampleRes);
			DiffuseDenoisedFBO2.SetSize(PADDED_WIDTH * DiffuseIndirectSuperSampleRes, PADDED_HEIGHT * DiffuseIndirectSuperSampleRes);
			VarianceFBO.SetSize(PADDED_WIDTH * DiffuseIndirectSuperSampleRes, PADDED_HEIGHT * DiffuseIndirectSuperSampleRes);


			if (TAA)
			{
				TAAFBO1.SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
				TAAFBO2.SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			}


			DownsampledFBO.SetSize(PADDED_WIDTH * 0.125f, PADDED_HEIGHT * 0.125f);
			BloomFBO.SetSize(PADDED_WIDTH * BloomQuality, PADDED_HEIGHT * BloomQuality);
			DiffractionSpikes.SetSize(PADDED_WIDTH * DiffractionResolution, PADDED_HEIGHT * DiffractionResolution);
			DiffractionSpikesDenoised.SetSize(PADDED_WIDTH * DiffractionResolution, PADDED_HEIGHT * DiffractionResolution);
			BloomFBOAlternate.SetSize(PADDED_WIDTH * BloomQuality, PADDED_HEIGHT * BloomQuality);
			BloomCombined.SetSize(PADDED_WIDTH * 0.5f, PADDED_HEIGHT * 0.5f);
			
			PostProcessingFBO.SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			DOFFBO.SetSize(TRUE_PADDED_WIDTH * DOFResolution, TRUE_PADDED_HEIGHT * DOFResolution);
			DOFFBOInput.SetSize(TRUE_PADDED_WIDTH * DOFResolution, TRUE_PADDED_HEIGHT * DOFResolution);
			TonemappedFBO.SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			FXAA_FBO.SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);

			//PostProcessingFBO_History[0].SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			//PostProcessingFBO_History[1].SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			//PostProcessingFBO_History[2].SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			//PostProcessingFBO_History[3].SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);
			//AntiFlickerFBO.SetSize(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);

			ColoredFBO.SetDimensions(TRUE_PADDED_WIDTH, TRUE_PADDED_HEIGHT);

			ShadowRawTrace.SetSize(PADDED_WIDTH * ShadowTraceResolution, PADDED_HEIGHT * ShadowTraceResolution);
			ShadowTemporalFBO_1.SetSize(PADDED_WIDTH * ShadowSupersampleRes, PADDED_HEIGHT * ShadowSupersampleRes);
			ShadowTemporalFBO_2.SetSize(PADDED_WIDTH * ShadowSupersampleRes, PADDED_HEIGHT * ShadowSupersampleRes);
			ShadowFiltered.SetSize(PADDED_WIDTH * ShadowSupersampleRes, PADDED_HEIGHT * ShadowSupersampleRes);
			ShadowSSS.SetSize(PADDED_WIDTH * ShadowTraceResolution, PADDED_HEIGHT * ShadowTraceResolution);
			ShadowSSS2.SetSize(PADDED_WIDTH * ShadowTraceResolution, PADDED_HEIGHT * ShadowTraceResolution);

			DenoiseSunShadows = DenoiseSunShadows && SoftShadows;

			ReflectionTraceFBO_1.SetSize(PADDED_WIDTH * ReflectionTraceResolution, PADDED_HEIGHT * ReflectionTraceResolution);
			ReflectionTraceFBO_2.SetSize(PADDED_WIDTH * ReflectionTraceResolution, PADDED_HEIGHT * ReflectionTraceResolution);
			ReflectionHitDataDenoised.SetSize(PADDED_WIDTH * ReflectionTraceResolution, PADDED_HEIGHT * ReflectionTraceResolution);
			
			ReflectionSuperSampleResolution = glm::max(ReflectionTraceResolution, ReflectionSuperSampleResolution);
			ReflectionTemporalFBO_1.SetSize(PADDED_WIDTH * ReflectionSuperSampleResolution, PADDED_HEIGHT * ReflectionSuperSampleResolution);
			ReflectionTemporalFBO_2.SetSize(PADDED_WIDTH * ReflectionSuperSampleResolution, PADDED_HEIGHT * ReflectionSuperSampleResolution);
			ReflectionDenoised_1.SetSize(PADDED_WIDTH * ReflectionSuperSampleResolution, PADDED_HEIGHT * ReflectionSuperSampleResolution);
			ReflectionDenoised_2.SetSize(PADDED_WIDTH * ReflectionSuperSampleResolution, PADDED_HEIGHT * ReflectionSuperSampleResolution);

			if (GodRays)
			{
				VolumetricFBO.SetSize(PADDED_WIDTH * VolumetricResolution, PADDED_HEIGHT * VolumetricResolution);
				BlurredVolumetricFBO.SetSize(PADDED_WIDTH * VolumetricResolution, PADDED_HEIGHT * VolumetricResolution);
			}

			if (SSAO)
			{
				SSAOFBO.SetSize(PADDED_WIDTH * SSAOResolution, PADDED_HEIGHT * SSAOResolution);
				SSAOBlurred.SetSize(PADDED_WIDTH * SSAOResolution, PADDED_HEIGHT * SSAOResolution);
			}

			if (RTAO)
			{
				float RTAO_Res2 = glm::max(RTAOResolution, 0.5f);
				RTAO_FBO.SetSize(PADDED_WIDTH * RTAOResolution, PADDED_HEIGHT * RTAOResolution);
				RTAO_TemporalFBO_1.SetSize(PADDED_WIDTH * RTAO_Res2, PADDED_HEIGHT * RTAO_Res2);
				RTAO_TemporalFBO_2.SetSize(PADDED_WIDTH * RTAO_Res2, PADDED_HEIGHT * RTAO_Res2);
			}
		}

		GLClasses::Framebuffer& TAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO1 : TAAFBO2;
		GLClasses::Framebuffer& PrevTAAFBO = (app.GetCurrentFrame() % 2 == 0) ? TAAFBO2 : TAAFBO1;
		GLClasses::Framebuffer& DiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO1 : DiffuseTemporalFBO2;
		GLClasses::Framebuffer& PrevDiffuseTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? DiffuseTemporalFBO2 : DiffuseTemporalFBO1;
		GLClasses::Framebuffer& ReflectionTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTemporalFBO_1 : ReflectionTemporalFBO_2;
		GLClasses::Framebuffer& PrevReflectionTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTemporalFBO_2 : ReflectionTemporalFBO_1;
		GLClasses::Framebuffer& RTAOTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? RTAO_TemporalFBO_1 : RTAO_TemporalFBO_2;
		GLClasses::Framebuffer& PrevRTAOTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? RTAO_TemporalFBO_2 : RTAO_TemporalFBO_1;
		GLClasses::Framebuffer& ShadowTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ShadowTemporalFBO_1 : ShadowTemporalFBO_2;
		GLClasses::Framebuffer& PrevShadowTemporalFBO = (app.GetCurrentFrame() % 2 == 0) ? ShadowTemporalFBO_2 : ShadowTemporalFBO_1;
		
		GLClasses::Framebuffer& ReflectionTraceFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTraceFBO_1 : ReflectionTraceFBO_2;
		GLClasses::Framebuffer& PrevReflectionTraceFBO = (app.GetCurrentFrame() % 2 == 0) ? ReflectionTraceFBO_2 : ReflectionTraceFBO_1;


		/////// history for anti flicker // unused.
		int framemod4 = app.GetCurrentFrame() % 4; // not true mathematical modulo, not needed since this wont be negative!
		//GLClasses::Framebuffer& PostProcessingFBO = PostProcessingFBO_History[framemod4];
		//GLClasses::Framebuffer& PostProcessingHistory0 = PostProcessingFBO_History[(framemod4+1)%4];
		//GLClasses::Framebuffer& PostProcessingHistory1 = PostProcessingFBO_History[(framemod4+2)%4];
		//GLClasses::Framebuffer& PostProcessingHistory2 = PostProcessingFBO_History[(framemod4+3)%4];



		if (glfwGetKey(app.GetWindow(), GLFW_KEY_F2) == GLFW_PRESS)
		{
			system("@cls");
			ShaderManager::RecompileShaders(); 
			world->m_ParticleEmitter.Recompile();
			Clouds::CloudRenderer::RecompileShaders();
			BloomRenderer::RecompileShaders();
			AtmosphereRenderer.Recompile();
			VoxelRT::Logger::Log("Recompiled!");
		}

		if (UpdatePlayerCollision && app.GetCurrentFrame() > 6) {
			bool TabPressed = glfwGetKey(app.GetWindow(), GLFW_KEY_TAB) == 1;
			MainPlayer.OnUpdate(app.GetWindow(), world, DeltaTime * 6.9f, (int)app.GetCurrentFrame(), DeltaSum, TabPressed);
		}

		MainPlayer.ClampVelocity();

		app.OnUpdate();

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
		float sun_visibility = glm::clamp(glm::dot(glm::normalize(SunDirection), glm::vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0f;

		// Projection update ->
		bool CloudProjectionTimeChangeUpdate = (CloudProjectionSunDir != StrongerLightDirection);

		if (app.GetCurrentFrame() %  12 == 0 || app.GetCurrentFrame() < 8 || CloudProjectionTimeChangeUpdate)
		{
			AtmosphereRenderer.RenderAtmosphere(SkymapMain, glm::normalize(SunDirection), 30, 4);
		}

		if (app.GetCurrentFrame() % 2 == 0 || CloudProjectionTimeChangeUpdate) {
			glm::mat4 CuberotationMap = glm::rotate(glm::mat4(1.0f), (float)glm::radians(glfwGetTime() * 2.0f), glm::vec3(0.3f, 0.1f, 1.0f));
			AtmosphereRenderer.DownsampleAtmosphere(SkymapSecondary, SkymapMain, CuberotationMap, NightSkyMapLowRes.GetID(), sun_visibility, NebulaGIColor ? NebulaGIStrength * 0.75f : 0.00f);
		}
		
		if (app.GetCurrentFrame() % 3 == 0 || CloudProjectionTimeChangeUpdate) {
			glm::mat4 CuberotationMap = glm::rotate(glm::mat4(1.0f), (float)glm::radians(glfwGetTime() * 2.0f), glm::vec3(0.3f, 0.1f, 1.0f));
			AtmosphereRenderer.DownsampleAtmosphere(SkymapSecondary_2, SkymapMain, CuberotationMap, NightSkyMapLowRes.GetID(), sun_visibility, NebulaGIColor ? NebulaGIStrength * 0.75f : 0.00f);
		}

		PreviousSunTick = SunTick;

		const int CloudProjectionUpdateRate = 359; 

		if ((app.GetCurrentFrame() == 4) || (app.GetCurrentFrame() == 12) || CloudProjectionTimeChangeUpdate
			|| (app.GetCurrentFrame() % CloudProjectionUpdateRate == 0) || (app.GetCurrentFrame() % (CloudProjectionUpdateRate+1) == 0) || (app.GetCurrentFrame() % (CloudProjectionUpdateRate+2) == 0)) {

			bool CutdownSteps = app.GetCurrentFrame() % CloudProjectionUpdateRate == 0 || app.GetCurrentFrame() % (CloudProjectionUpdateRate + 1) == 0 || app.GetCurrentFrame() % (CloudProjectionUpdateRate + 2) == 0;



			// Cube views ->
			glm::vec3 VirtualCenter = glm::vec3(MainPlayer.m_Position.x, 150.0, MainPlayer.m_Position.z);
			auto Matrices = GetMatrices(VirtualCenter);
			glm::mat4 VirtualProjection = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, 1000.0f);

			int StartIndex = 0;

			if (CutdownSteps) {
				if (app.GetCurrentFrame() % CloudProjectionUpdateRate == 0) {
					StartIndex = 0;
				}

				else if (app.GetCurrentFrame() % (CloudProjectionUpdateRate + 1) == 0) {
					StartIndex = 2;
				}

				else if (app.GetCurrentFrame() % (CloudProjectionUpdateRate + 2) == 0) {
					StartIndex = 4;
				}
			}

			glm::vec3 CelestialDirections[2] = { SunDirection, MoonDirection };

			// Render clouds for the entire cubemap ->
			for (int i = StartIndex; i < 6; i++) {

				glm::mat4 ViewMatrixCube = Matrices[i];

				GLuint CloudFramebuffer = Clouds::CloudRenderer::Update(VirtualProjection, ViewMatrixCube, VirtualProjection, ViewMatrixCube, VirtualCenter,
					VirtualCenter, VAO, StrongerLightDirection, BluenoiseTexture.GetTextureID(),
					PADDED_WIDTH, PADDED_HEIGHT, app.GetCurrentFrame(), SkymapSecondary.GetTexture(), InitialTraceFBO->GetTexture(0), VirtualCenter, InitialTraceFBOPrev->GetTexture(0),
					CloudModifiers, ClampCloudTemporal, glm::vec3(CloudDetailScale, CloudDetailWeightEnabled ? 1.0f : 0.0f, CloudErosionWeightExponent),
					CloudTimeScale, CurlNoiseOffset, CirrusStrength, CirrusScale, app.GetCurrentFrame() == 4 ? glm::ivec3(64, 8, 4) : (CutdownSteps ? glm::ivec3(32, 6, 3) : glm::ivec3(32, 6, 3)), CloudCheckerStepCount, sun_visibility, CloudDetailFBMPower,
					CloudLODLighting, CloudForceSupersample, CloudForceSupersampleRes, CloudSpatialUpscale, CloudAmbientDensityMultiplier, CloudProjection.GetTexture(0), true, CloudThiccness, CloudDetailScale, true, CelestialDirections, SkymapMain.GetTexture(), glm::vec2(SunStrengthModifier, MoonStrengthModifier));


				glBindFramebuffer(GL_FRAMEBUFFER, CloudFramebuffer);
				glClear(GL_COLOR_BUFFER_BIT);
			}


			std::cout << "\n\n---CLOUD PROJECTION UPDATE---\n\n";
			std::cout << "\nCLOUD PROJECTION DETAILS : \n";
			std::cout << "\nTIME UPDATE : " << CutdownSteps << "\tSTART INDEX : " << StartIndex << "\tSUN_VISIB : " << sun_visibility;
			std::cout << "\n";
		}

		CloudProjectionSunDir = StrongerLightDirection;








		// GBuffer ->

		bool UpdateGBufferThisFrame = PreviousView != CurrentView || app.GetCurrentFrame() % 10 == 0 ||
			ModifiedWorld || JitterSceneForTAA;

		if (UpdateGBufferThisFrame)
		{
			if (HighlightFocusedBlock) {
				FocusedOnBlock = world->RaycastDetect(MainCamera.GetPosition(), MainCamera.GetFront());
			}

			// Swap the initial trace framebuffers
			InitialTraceFBO = InitialTraceFBO == &InitialTraceFBO_1 ? &InitialTraceFBO_2 : &InitialTraceFBO_1;
			InitialTraceFBOPrev = InitialTraceFBO == &InitialTraceFBO_1 ? &InitialTraceFBO_2 : &InitialTraceFBO_1;

			InitialTraceFBO->Bind();
			glDisable(GL_CULL_FACE);
			glDisable(GL_DEPTH_TEST);

			InitialTraceShader.Use();

			glm::mat4 JitterMatrix = GetTAAJitterMatrix(app.GetCurrentFrame(), glm::vec2(floor(app.GetWidth() * InitialTraceResolution), floor(app.GetHeight() * InitialTraceResolution)));
			glm::vec2 TAAJitter = GetTAAJitter(app.GetCurrentFrame(), glm::vec2(floor(PADDED_WIDTH * InitialTraceResolution), floor(PADDED_HEIGHT * InitialTraceResolution)));


			InitialTraceShader.SetMatrix4("u_InverseView", inv_view);
			InitialTraceShader.SetMatrix4("u_InverseProjection",  glm::inverse(MainCamera.GetProjectionMatrix()));
			InitialTraceShader.SetInteger("u_VoxelDataTexture", 0);
			InitialTraceShader.SetInteger("u_AlbedoTextures", 1);
			InitialTraceShader.SetInteger("u_RenderDistance", RenderDistance);
			InitialTraceShader.SetInteger("u_DistanceFieldTexture", 2);
			InitialTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			InitialTraceShader.SetInteger("u_VertCurrentFrame", app.GetCurrentFrame());
			InitialTraceShader.SetVector2f("u_Dimensions", glm::vec2(InitialTraceFBO->GetWidth(), InitialTraceFBO->GetHeight()));
			InitialTraceShader.SetVector2f("u_VertDimensions", glm::vec2(PADDED_WIDTH, PADDED_HEIGHT));
			InitialTraceShader.SetVector2f("u_CurrentTAAJitter", glm::vec2(TAAJitter));
			InitialTraceShader.SetVector3f("u_PlayerPosition", MainCamera.GetPosition());
			InitialTraceShader.SetFloat("u_FOV", MainCamera.GetFov());
			InitialTraceShader.SetFloat("u_TanFOV", glm::tan(MainCamera.GetFov()));
			InitialTraceShader.SetFloat("u_Time", glfwGetTime());
			InitialTraceShader.SetBool("u_ShouldAlphaTest", ShouldAlphaTest);
			InitialTraceShader.SetBool("u_JitterSceneForTAA", JitterSceneForTAA);
			
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

			BlockDataStorageBuffer.Bind(0);

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			InitialTraceFBO->Unbind();

			if (DOWNSAMPLE_GBUFFERS) {
				// Downsampling ->

				// Half res ->
				GBuffer4xDownsample.Use();
				HalfResGBuffer.Bind();
				GBuffer4xDownsample.SetInteger("u_DepthFullRes", 0);
				GBuffer4xDownsample.SetInteger("u_NormalsFullRes", 1);
				GBuffer4xDownsample.SetInteger("u_BlockIDs", 2);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();

				// Quarter res ->
				QuarterResGBuffer.Bind();
				GBuffer4xDownsample.Use();
				GBuffer4xDownsample.SetInteger("u_DepthFullRes", 0);
				GBuffer4xDownsample.SetInteger("u_NormalsFullRes", 1);
				GBuffer4xDownsample.SetInteger("u_BlockIDs", 2);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, HalfResGBuffer.GetTexture(0));

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, HalfResGBuffer.GetTexture(1));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, HalfResGBuffer.GetTexture(2));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}
			
		}



		{
			// Generate gbuffer ->

			GenerateGBuffer.Use();
			GeneratedGBuffer.Bind();

			GenerateGBuffer.SetInteger("u_NonLinearDepth", 0);
			GenerateGBuffer.SetInteger("u_Normals", 1);
			GenerateGBuffer.SetInteger("u_BlockIDs", 2);
			GenerateGBuffer.SetInteger("u_BlockAlbedos", 3);
			GenerateGBuffer.SetInteger("u_BlockNormals", 4);
			GenerateGBuffer.SetInteger("u_BlockPBR", 5);
			GenerateGBuffer.SetInteger("u_BlockEmissive", 6);

			GenerateGBuffer.SetInteger("u_LavaTextures[0]", 7);
			GenerateGBuffer.SetInteger("u_LavaTextures[1]", 8);

			GenerateGBuffer.SetMatrix4("u_InverseView", inv_view);
			GenerateGBuffer.SetMatrix4("u_InverseProjection", glm::inverse(MainCamera.GetProjectionMatrix()));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[0]", VoxelRT::BlockDatabase::GetBlockID("Grass"));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[1]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[2]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[3]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[4]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[5]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[6]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[7]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[8]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			GenerateGBuffer.SetInteger("u_GrassBlockProps[9]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			GenerateGBuffer.SetFloat("u_Time", glfwGetTime());
			GenerateGBuffer.SetFloat("uTime", glfwGetTime());
			GenerateGBuffer.SetInteger("u_LavaBlockID", BlockDatabase::GetBlockID("Lava"));
			GenerateGBuffer.SetBool("u_POM", POM);
			GenerateGBuffer.SetFloat("u_POMHeight", POMHeight);
			GenerateGBuffer.SetBool("u_HighQualityPOM", HighQualityPOM);
			GenerateGBuffer.SetBool("u_DitherPOM", DitherPOM);
			GenerateGBuffer.SetBool("u_UpdateGBufferThisFrame", UpdateGBufferThisFrame);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

			glActiveTexture(GL_TEXTURE6);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

			glActiveTexture(GL_TEXTURE7);
			glBindTexture(GL_TEXTURE_3D, LavaAlbedo.m_ID);

			glActiveTexture(GL_TEXTURE8);
			glBindTexture(GL_TEXTURE_3D, LavaNormals.m_ID);

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}



		if (DOF) {
			WriteCenterDepth.Use();

			WriteCenterDepth.SetInteger("u_DepthTexture", 0);
			WriteCenterDepth.SetVector2f("u_CameraPlanes", glm::vec2(MainCamera.GetNearPlane(), MainCamera.GetFarPlane()));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

			glBindBuffer(GL_SHADER_STORAGE_BUFFER, DOFSSBO);
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, DOFSSBO);

			glDispatchCompute(1, 1, 1);

			glUseProgram(0);

			float DownloadedCenterDepth = 0.0f;

			glBindBuffer(GL_SHADER_STORAGE_BUFFER, DOFSSBO);
			glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(float), &DownloadedCenterDepth);;
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

			CenterDepthSmooth = glm::mix(DownloadedCenterDepth, CenterDepthSmooth, DOFTemporalDepthBlend);

			if (app.GetCurrentFrame() % 16 == 0)
				std::cout << "\n\nDOF Temporal Center Depth : " << CenterDepthSmooth << "\n";
		}




		// Diffuse tracing

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		DiffuseRawTraceFBO.Bind();

		DiffuseTraceShader.Use();
		
		DiffuseTraceShader.SetInteger("u_VoxelData", 0);
		DiffuseTraceShader.SetInteger("u_PositionTexture", 1);
		DiffuseTraceShader.SetInteger("u_NormalTexture", 2);
		DiffuseTraceShader.SetInteger("u_Skymap", 3);
		DiffuseTraceShader.SetInteger("u_BlockNormalTextures", 4);
		DiffuseTraceShader.SetInteger("u_BlockAlbedoTextures", 6);
		DiffuseTraceShader.SetInteger("u_BlueNoiseTextures", 7);
		DiffuseTraceShader.SetInteger("u_BlockPBRTextures", 8);
		DiffuseTraceShader.SetInteger("u_BlockEmissiveTextures", 11);
		DiffuseTraceShader.SetInteger("u_DistanceFieldTexture", 13);
		DiffuseTraceShader.SetInteger("u_DiffuseTraceLength", DiffuseTraceLength);
		DiffuseTraceShader.SetVector2f("u_Halton", glm::vec2(GetTAAJitterSecondary(app.GetCurrentFrame())));

		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
		DiffuseTraceShader.SetMatrix4("u_InverseView", inv_view);
		DiffuseTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
		DiffuseTraceShader.SetVector2f("u_Dimensions", glm::vec2(DiffuseRawTraceFBO.GetWidth(), DiffuseRawTraceFBO.GetHeight()));
		DiffuseTraceShader.SetFloat("u_Time", glfwGetTime() * 1.2f);
		DiffuseTraceShader.SetFloat("u_DiffuseLightIntensity", DiffuseLightIntensity);

		
		if (app.GetCurrentFrame() % 100 == 0) {
			std::cout << "\nSUN VISIBILITY : " << sun_visibility;
		}

		DiffuseTraceShader.SetFloat("u_SunVisibility", sun_visibility);
		DiffuseTraceShader.SetFloat("u_GISunStrength", GISunStrength);
		DiffuseTraceShader.SetFloat("u_GISkyStrength", GISkyStrength);
		DiffuseTraceShader.SetFloat("u_SunVisibility", sun_visibility);
		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
		DiffuseTraceShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
		DiffuseTraceShader.SetVector3f("u_SunDirection", glm::normalize(SunDirection));
		DiffuseTraceShader.SetVector3f("u_MoonDirection", glm::normalize(MoonDirection));

		DiffuseTraceShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
		DiffuseTraceShader.SetMatrix4("u_ShadowView", ShadowView);
		DiffuseTraceShader.SetInteger("u_ShadowMap", 9);
		DiffuseTraceShader.SetInteger("u_BlueNoiseTexture", 10);
		DiffuseTraceShader.SetInteger("u_SPP", DiffuseSPP);
		DiffuseTraceShader.SetInteger("u_CheckerSPP", (DiffuseSPP + DiffuseSPP % 2) / 2); 
		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());

		DiffuseTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
		DiffuseTraceShader.SetInteger("u_CurrentFrameMod512", glm::clamp((int)app.GetCurrentFrame() % 512, 0, 524));

		DiffuseTraceShader.SetInteger("u_CurrentFrameMod128", (int)app.GetCurrentFrame() % 128);
		DiffuseTraceShader.SetBool("u_UseBlueNoise", USE_BLUE_NOISE_FOR_TRACING);
		DiffuseTraceShader.SetBool("CHECKERBOARD_SPP", CHECKERBOARD_SPP);
		DiffuseTraceShader.SetBool("u_APPLY_PLAYER_SHADOW", APPLY_PLAYER_SHADOW_FOR_GI);
		DiffuseTraceShader.SetBool("u_DirectSampling", DiffuseDirectLightSampling);
		DiffuseTraceShader.SetBool("u_UseDirectSampling", DiffuseDirectLightSampling);
		DiffuseTraceShader.SetBool("u_UseMIS", DiffuseDirectLightSampling);
		DiffuseTraceShader.SetBool("u_Supersample", abs(DiffuseIndirectSuperSampleRes - DiffuseTraceResolution) > 0.05f);
		DiffuseTraceShader.SetMatrix4("u_VertInverseView", inv_view);
		DiffuseTraceShader.SetMatrix4("u_VertInverseProjection", inv_projection);
		DiffuseTraceShader.SetMatrix4("u_InverseView", inv_view);
		DiffuseTraceShader.SetMatrix4("u_InverseProjection", inv_projection);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapSecondary_2.GetTexture());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlueNoise.GetTextureArray());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, ShadowRawTrace.GetTexture());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

		glActiveTexture(GL_TEXTURE13);
		glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

		BlockDataStorageBuffer.Bind(0);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, BlueNoise_SSBO.m_SSBO);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, world->LightChunkDataSSBO);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, world->LightChunkOffsetSSBO);

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		DiffuseRawTraceFBO.Unbind();


		if (USE_SVGF)
		{
				// Initial 3x3 atrous filter 

				if (PreTemporalSpatialPass) {

					if (app.GetCurrentFrame() % 32 == 0) { std::cout << "\n! pre spatial pass !\n"; }

					SpatialInitial.Use();
					DiffusePreTemporal_SpatialFBO.Bind();


					SpatialInitial.SetInteger("u_SH", 0);
					SpatialInitial.SetInteger("u_CoCg", 1);
					SpatialInitial.SetInteger("u_Utility", 2);
					SpatialInitial.SetInteger("u_AO", 3);
					SpatialInitial.SetInteger("u_PositionTexture", 4);
					SpatialInitial.SetInteger("u_NormalTexture", 5);


					SpatialInitial.SetVector2f("u_Dimensions", glm::vec2(DiffusePreTemporal_SpatialFBO.GetWidth(), DiffusePreTemporal_SpatialFBO.GetHeight()));
					SpatialInitial.SetMatrix4("u_VertInverseView", inv_view);
					SpatialInitial.SetMatrix4("u_VertInverseProjection", inv_projection);
					SpatialInitial.SetMatrix4("u_InverseView", inv_view);
					SpatialInitial.SetMatrix4("u_InverseProjection", inv_projection);
					SpatialInitial.SetFloat("u_Time", glfwGetTime());
					SpatialInitial.SetFloat("u_DeltaTime", DeltaTime);

					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(0));
					glActiveTexture(GL_TEXTURE1);
					glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(1));
					glActiveTexture(GL_TEXTURE2);
					glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(2));
					glActiveTexture(GL_TEXTURE3);
					glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(3));

					// Gbuffer 
					glActiveTexture(GL_TEXTURE4);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));
					glActiveTexture(GL_TEXTURE5);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

					VAO.Bind();
					glDrawArrays(GL_TRIANGLES, 0, 6);
					VAO.Unbind();
				}



				// Temporal filter 

				DiffuseTemporalFBO.Bind();
				SVGF_Temporal.Use();

				SVGF_Temporal.SetInteger("u_CurrentCoCg", 0);
				SVGF_Temporal.SetInteger("u_CurrentPositionTexture", 1);
				SVGF_Temporal.SetInteger("u_PrevCoCg", 2);
				SVGF_Temporal.SetInteger("u_PreviousPositionTexture", 3);


				SVGF_Temporal.SetInteger("u_CurrentSH", 4);
				SVGF_Temporal.SetInteger("u_PreviousSH", 5);

				SVGF_Temporal.SetInteger("u_CurrentNormalTexture", 6);
				SVGF_Temporal.SetInteger("u_PreviousNormalTexture", 7);
				SVGF_Temporal.SetInteger("u_PreviousUtility", 8);
				SVGF_Temporal.SetInteger("u_NoisyLuminosity", 9);

				SVGF_Temporal.SetInteger("u_CurrentBlockIDTexture", 10);
				SVGF_Temporal.SetInteger("u_PrevBlockIDTexture", 11);


				//u_CurrentAO;
				//u_PreviousAO;
				SVGF_Temporal.SetInteger("u_CurrentAO", 12);
				SVGF_Temporal.SetInteger("u_PreviousAO", 13);


				SVGF_Temporal.SetBool("u_DiffuseTemporal", true);
				SVGF_Temporal.SetBool("u_ShadowTemporal", false);


				SVGF_Temporal.SetMatrix4("u_Projection", CurrentProjection);
				SVGF_Temporal.SetMatrix4("u_View", CurrentView);
				SVGF_Temporal.SetMatrix4("u_PrevProjection", PreviousProjection);
				SVGF_Temporal.SetMatrix4("u_PrevView", PreviousView);

				SVGF_Temporal.SetFloat("u_MinimumMix", 0.0f);
				SVGF_Temporal.SetFloat("u_MaximumMix", 0.96f);
				SVGF_Temporal.SetInteger("u_TemporalQuality", 0); // No clamping!
				SVGF_Temporal.SetBool("u_ReflectionTemporal", false);
				SVGF_Temporal.SetBool("u_BeUseful", DO_SVGF_TEMPORAL);
				SVGF_Temporal.SetFloat("u_ClampBias", 0.025f);
				SVGF_Temporal.SetVector3f("u_PrevCameraPos", PreviousPosition);
				SVGF_Temporal.SetVector3f("u_CurrentCameraPos", MainCamera.GetPosition());
				SVGF_Temporal.SetVector2f("u_Dimensions", glm::vec2(DiffuseTemporalFBO.GetWidth(), DiffuseTemporalFBO.GetHeight()));


				SVGF_Temporal.SetMatrix4("u_VertInverseView", inv_view);
				SVGF_Temporal.SetMatrix4("u_VertInverseProjection", inv_projection);
				SVGF_Temporal.SetMatrix4("u_InverseView", inv_view);
				SVGF_Temporal.SetMatrix4("u_InverseProjection", inv_projection);
				SVGF_Temporal.SetMatrix4("u_PrevInverseProjection", glm::inverse(PreviousProjection));
				SVGF_Temporal.SetMatrix4("u_PrevInverseView", glm::inverse(PreviousView));

				SVGF_Temporal.SetFloat("u_Time", glfwGetTime());
				SVGF_Temporal.SetFloat("u_DeltaTime", DeltaTime);

				//uniform mat4 u_PrevInverseProjection;
				//uniform mat4 u_PrevInverseView;
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, PreTemporalSpatialPass ? DiffusePreTemporal_SpatialFBO.GetTexture(1) : DiffuseRawTraceFBO.GetTexture(1));
				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture(1));
				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(0));

				glActiveTexture(GL_TEXTURE4);
				glBindTexture(GL_TEXTURE_2D, PreTemporalSpatialPass ? DiffusePreTemporal_SpatialFBO.GetTexture(0) : DiffuseRawTraceFBO.GetTexture(0));
				glActiveTexture(GL_TEXTURE5);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture(0));

				glActiveTexture(GL_TEXTURE6);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));
				glActiveTexture(GL_TEXTURE7);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(1));

				glActiveTexture(GL_TEXTURE8);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture(2));

				glActiveTexture(GL_TEXTURE9);
				glBindTexture(GL_TEXTURE_2D, PreTemporalSpatialPass ? DiffusePreTemporal_SpatialFBO.GetTexture(2) : DiffuseRawTraceFBO.GetTexture(2));

				glActiveTexture(GL_TEXTURE10);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

				glActiveTexture(GL_TEXTURE11);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(2));

				glActiveTexture(GL_TEXTURE12);
				glBindTexture(GL_TEXTURE_2D, PreTemporalSpatialPass ? DiffusePreTemporal_SpatialFBO.GetTexture(3) : DiffuseRawTraceFBO.GetTexture(3));
				glActiveTexture(GL_TEXTURE13);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture(3));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();

				DiffuseTemporalFBO.Unbind();

				// Do a variance estimation pass

				VarianceFBO.Bind();
				SVGF_Variance.Use();

				SVGF_Variance.SetInteger("u_PositionTexture", 0);
				SVGF_Variance.SetInteger("u_NormalTexture", 1);
				SVGF_Variance.SetInteger("u_SH", 2);
				SVGF_Variance.SetInteger("u_CoCg", 3);
				SVGF_Variance.SetInteger("u_Utility", 4);
				SVGF_Variance.SetMatrix4("u_InverseView", inv_view);
				SVGF_Variance.SetMatrix4("u_InverseProjection", inv_projection);
				SVGF_Variance.SetMatrix4("u_VertInverseView", inv_view);
				SVGF_Variance.SetMatrix4("u_VertInverseProjection", inv_projection);
				SVGF_Variance.SetBool("DO_SPATIAL", DO_VARIANCE_SPATIAL);
				SVGF_Variance.SetBool("AGGRESSIVE_DISOCCLUSION_HANDLING", AGGRESSIVE_DISOCCLUSION_HANDLING);
				SVGF_Variance.SetBool("u_LargeKernel", SVGF_LARGE_KERNEL);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture());

				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture(1));

				glActiveTexture(GL_TEXTURE4);
				glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture(2));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();

				VarianceFBO.Unbind();

				// Spatial filter :

				int StepSizes[5];


				if (WiderSVGF)
				{
					StepSizes[0] = 32;
					StepSizes[1] = 16;
					StepSizes[2] = 8;
					StepSizes[3] = 4;
					StepSizes[4] = 2;
				}
				
				else
				{
					StepSizes[0] = 12;
					StepSizes[1] = 8;
					StepSizes[2] = 6;
					StepSizes[3] = 4;
					StepSizes[4] = 2;
				}

				for (int i = 0; i < 5; i++)
				{
					// 1 2 1 2 1
					auto& CurrentDenoiseFBO = (i % 2 == 0) ? DiffuseDenoiseFBO : DiffuseDenoisedFBO2;
					auto& PrevDenoiseFBO = (i == 0) ? VarianceFBO :
						(i % 2 == 0) ? DiffuseDenoisedFBO2 : DiffuseDenoiseFBO;

					GLuint VarianceTexture = 0;
					GLuint AOTexture = 0;

					if (i == 0)
					{
						VarianceTexture = VarianceFBO.GetTexture(2);
					}

					else {

						if (i % 2 == 0)
						{
							VarianceTexture = DiffuseDenoisedFBO2.GetTexture(2);
						}

						else {
							VarianceTexture = DiffuseDenoiseFBO.GetTexture(2);
						}
					}

					// ao texture

					if (i == 0)
					{
						AOTexture = DiffuseTemporalFBO.GetTexture(3);
					}

					else {

						if (i % 2 == 0)
						{
							AOTexture = DiffuseDenoisedFBO2.GetTexture(3);
						}

						else {
							AOTexture = DiffuseDenoiseFBO.GetTexture(3);
						}
					}







					CurrentDenoiseFBO.Bind();
					SVGF_Spatial.Use();

					// textures :
					SVGF_Spatial.SetInteger("u_SH", 0);
					SVGF_Spatial.SetInteger("u_PositionTexture", 1);
					SVGF_Spatial.SetInteger("u_NormalTexture", 2);
					SVGF_Spatial.SetInteger("u_BlockIDTexture", 3);
					SVGF_Spatial.SetInteger("u_VarianceTexture", 4);
					SVGF_Spatial.SetInteger("u_CoCg", 5);
					SVGF_Spatial.SetInteger("u_Utility", 6);
					SVGF_Spatial.SetInteger("u_AO", 8);
					SVGF_Spatial.SetInteger("u_TemporalMoment", 9);
					SVGF_Spatial.SetInteger("u_CurrentPass", i);
					SVGF_Spatial.SetBool("u_LargeKernel", SVGF_LARGE_KERNEL);

					SVGF_Spatial.SetInteger("u_Step", StepSizes[i]);
					SVGF_Spatial.SetVector2f("u_Dimensions", glm::vec2(CurrentDenoiseFBO.GetWidth(), CurrentDenoiseFBO.GetHeight()));
					SVGF_Spatial.SetMatrix4("u_VertInverseView", inv_view);
					SVGF_Spatial.SetMatrix4("u_VertInverseProjection", inv_projection);
					SVGF_Spatial.SetMatrix4("u_InverseView", inv_view);
					SVGF_Spatial.SetMatrix4("u_InverseProjection", inv_projection);
					SVGF_Spatial.SetBool("u_ShouldDetailWeight", !(i >= 3));
					SVGF_Spatial.SetBool("DO_SPATIAL", DO_SVGF_SPATIAL);
					SVGF_Spatial.SetBool("AGGRESSIVE_DISOCCLUSION_HANDLING", AGGRESSIVE_DISOCCLUSION_HANDLING);
					SVGF_Spatial.SetFloat("u_ColorPhiBias", ColorPhiBias);
					SVGF_Spatial.SetFloat("u_Time", glfwGetTime());
					SVGF_Spatial.SetFloat("u_DeltaTime", DeltaTime);
					SVGF_Spatial.SetFloat("u_ResolutionScale", DiffuseIndirectSuperSampleRes);

					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, PrevDenoiseFBO.GetTexture());

					glActiveTexture(GL_TEXTURE1);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

					glActiveTexture(GL_TEXTURE2);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

					glActiveTexture(GL_TEXTURE3);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

					glActiveTexture(GL_TEXTURE4);
					glBindTexture(GL_TEXTURE_2D, VarianceTexture);

					glActiveTexture(GL_TEXTURE5);
					glBindTexture(GL_TEXTURE_2D, PrevDenoiseFBO.GetTexture(1));

					glActiveTexture(GL_TEXTURE6);
					glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture(3));

					glActiveTexture(GL_TEXTURE8);
					glBindTexture(GL_TEXTURE_2D, AOTexture);

					glActiveTexture(GL_TEXTURE9);
					glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture(2));

					VAO.Bind();
					glDrawArrays(GL_TRIANGLES, 0, 6);
					VAO.Unbind();
				}

				glUseProgram(0);
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}


		// ATROUS WAVELET FILTER //
		else
		{
				// Temporal filter 

				DiffuseTemporalFBO.Bind();
				MainTemporalFilter.Use();

				MainTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
				MainTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
				MainTemporalFilter.SetInteger("u_PreviousColorTexture", 2);
				MainTemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);


				MainTemporalFilter.SetInteger("u_CurrentSH", 4);
				MainTemporalFilter.SetInteger("u_PreviousSH", 5);
				MainTemporalFilter.SetInteger("u_NormalTexture", 8);

				MainTemporalFilter.SetInteger("u_CurrentAO", 10);
				MainTemporalFilter.SetInteger("u_PreviousAO", 11);

				MainTemporalFilter.SetBool("u_DiffuseTemporal", true);
				MainTemporalFilter.SetBool("u_ShadowTemporal", false);


				MainTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
				MainTemporalFilter.SetMatrix4("u_View", CurrentView);
				MainTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
				MainTemporalFilter.SetMatrix4("u_PrevView", PreviousView);

				MainTemporalFilter.SetFloat("u_MinimumMix", 0.0f);
				MainTemporalFilter.SetFloat("u_MaximumMix", 0.96f);
				MainTemporalFilter.SetInteger("u_TemporalQuality", 0); // No clamping!
				MainTemporalFilter.SetBool("u_ReflectionTemporal", false);
				MainTemporalFilter.SetFloat("u_ClampBias", 0.025f);
				MainTemporalFilter.SetVector3f("u_PrevCameraPos", PreviousPosition);
				MainTemporalFilter.SetVector3f("u_CurrentCameraPos", MainCamera.GetPosition());


				MainTemporalFilter.SetMatrix4("u_VertInverseView", inv_view);
				MainTemporalFilter.SetMatrix4("u_VertInverseProjection", inv_projection);
				MainTemporalFilter.SetMatrix4("u_InverseView", inv_view);
				MainTemporalFilter.SetMatrix4("u_InverseProjection", inv_projection);
				MainTemporalFilter.SetMatrix4("u_PrevInverseProjection", glm::inverse(PreviousProjection));
				MainTemporalFilter.SetMatrix4("u_PrevInverseView", glm::inverse(PreviousView));

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture());

				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(0));

				glActiveTexture(GL_TEXTURE4);
				glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(1));

				glActiveTexture(GL_TEXTURE5);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture(1));

				glActiveTexture(GL_TEXTURE8);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

				glActiveTexture(GL_TEXTURE10);
				glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(3));
				glActiveTexture(GL_TEXTURE11);
				glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporalFBO.GetTexture(3));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();

				DiffuseTemporalFBO.Unbind();

				// Do a variance estimation pass

				VarianceFBO.Bind();
				VarianceEstimator.Use();
				VarianceEstimator.SetInteger("u_InputTexture", 0);
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, DiffuseTemporalFBO.GetTexture());
				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
				VarianceFBO.Unbind();

				//
				/// Do 5 atrous spatial filter pass with varying step sizes
				//

				int StepSizes[5] = { 12, 8, 4, 2, 1 };

				for (int i = 0; i < 5; i++)
				{
					GLuint AOTexture = 0;

					// 1 2 1 2 1
					auto& CurrentDenoiseFBO = (i % 2 == 0) ? DiffuseDenoiseFBO : DiffuseDenoisedFBO2;
					auto& PrevDenoiseFBO = (i == 0) ? DiffuseTemporalFBO :
						(i % 2 == 0) ? DiffuseDenoisedFBO2 : DiffuseDenoiseFBO;

					CurrentDenoiseFBO.Bind();
					AtrousSpatialFilter.Use();
					AtrousSpatialFilter.SetInteger("u_InputTexture", 0);
					AtrousSpatialFilter.SetInteger("u_PositionTexture", 1);
					AtrousSpatialFilter.SetInteger("u_NormalTexture", 2);
					AtrousSpatialFilter.SetInteger("u_BlockIDTexture", 3);
					AtrousSpatialFilter.SetInteger("u_VarianceTexture", 4);
					AtrousSpatialFilter.SetInteger("u_InputTexture2", 5);
					AtrousSpatialFilter.SetInteger("u_AO", 6);
					AtrousSpatialFilter.SetInteger("u_Step", StepSizes[i]);
					AtrousSpatialFilter.SetVector2f("u_Dimensions", glm::vec2(CurrentDenoiseFBO.GetWidth(), CurrentDenoiseFBO.GetHeight()));
					AtrousSpatialFilter.SetMatrix4("u_VertInverseView", inv_view);
					AtrousSpatialFilter.SetMatrix4("u_VertInverseProjection", inv_projection);
					AtrousSpatialFilter.SetMatrix4("u_InverseView", inv_view);
					AtrousSpatialFilter.SetMatrix4("u_InverseProjection", inv_projection);
					AtrousSpatialFilter.SetBool("u_ShouldDetailWeight", !(i >= 3));

					// ao texture

					if (i == 0)
					{
						AOTexture = DiffuseTemporalFBO.GetTexture(3);
					}

					else {

						if (i % 2 == 0)
						{
							AOTexture = DiffuseDenoisedFBO2.GetTexture(3);
						}

						else {
							AOTexture = DiffuseDenoiseFBO.GetTexture(3);
						}
					}

					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, PrevDenoiseFBO.GetTexture());

					glActiveTexture(GL_TEXTURE1);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

					glActiveTexture(GL_TEXTURE2);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

					glActiveTexture(GL_TEXTURE3);
					glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

					glActiveTexture(GL_TEXTURE4);
					glBindTexture(GL_TEXTURE_2D, VarianceFBO.GetTexture());

					glActiveTexture(GL_TEXTURE5);
					glBindTexture(GL_TEXTURE_2D, PrevDenoiseFBO.GetTexture(1));

					glActiveTexture(GL_TEXTURE6);
					glBindTexture(GL_TEXTURE_2D, AOTexture);

					VAO.Bind();
					glDrawArrays(GL_TRIANGLES, 0, 6);
					VAO.Unbind();
				}

				glUseProgram(0);
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}


		// ---- SHADOW TRACE ----

		GLClasses::Framebuffer& ShadowFBO = ShadowRawTrace;

		{
			ShadowFBO.Bind();
			ShadowTraceShader.Use();

			ShadowTraceShader.SetInteger("u_PositionTexture", 0);
			ShadowTraceShader.SetInteger("u_VoxelData", 1);
			ShadowTraceShader.SetInteger("u_AlbedoTextures", 2);
			ShadowTraceShader.SetInteger("u_NormalTexture", 3);
			ShadowTraceShader.SetInteger("u_DistanceFieldTexture", 5);
			ShadowTraceShader.SetInteger("u_BlueNoiseTexture", 6);

			ShadowTraceShader.SetVector3f("u_LightDirection", StrongerLightDirection);
			ShadowTraceShader.SetVector3f("u_PlayerPosition", MainCamera.GetPosition());
			ShadowTraceShader.SetVector2f("u_Dimensions", glm::vec2(ShadowFBO.GetWidth(), ShadowFBO.GetHeight()));
			ShadowTraceShader.SetBool("u_DoFullTrace", true);
			ShadowTraceShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
			ShadowTraceShader.SetMatrix4("u_ShadowView", ShadowView);

			ShadowTraceShader.SetMatrix4("u_VertInverseView", inv_view);
			ShadowTraceShader.SetMatrix4("u_VertInverseProjection", inv_projection);
			ShadowTraceShader.SetMatrix4("u_InverseView", inv_view);
			ShadowTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
			ShadowTraceShader.SetFloat("u_Time", glfwGetTime());
			ShadowTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());

			ShadowTraceShader.SetBool("u_ContactHardeningShadows", SoftShadows);
			ShadowTraceShader.SetBool("u_ShouldAlphaTest", ShouldAlphaTestShadows);
			ShadowTraceShader.SetFloat("u_FOV", MainCamera.GetFov());
			ShadowTraceShader.SetFloat("u_TanFOV", glm::tan(MainCamera.GetFov()));
			ShadowTraceShader.SetVector2f("u_Halton", glm::abs(ShadowSupersampleRes - ShadowTraceResolution) > 0.1f ? glm::vec2(GetTAAJitterSecondary(app.GetCurrentFrame())) : glm::vec2(0.));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE6);
			glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

			BlockDataStorageBuffer.Bind(0);

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ShadowFBO.Unbind();
			ShadowProjection = CurrentProjection;
			ShadowView = CurrentView;
		}

		if (SoftShadows)
		{
			ShadowTemporalFBO.Bind();

			ShadowTemporalFilter.Use();
			ShadowTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
			ShadowTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
			ShadowTemporalFilter.SetInteger("u_PreviousColorTexture", 2);
			ShadowTemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);
			ShadowTemporalFilter.SetInteger("u_NormalTexture", 8);
			ShadowTemporalFilter.SetInteger("u_ShadowTransversals", 17);
			ShadowTemporalFilter.SetInteger("u_DenoisedTransversals", 19);
			ShadowTemporalFilter.SetInteger("u_FrameCount", 22);

			ShadowTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
			ShadowTemporalFilter.SetMatrix4("u_View", CurrentView);
			ShadowTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
			ShadowTemporalFilter.SetMatrix4("u_PrevView", PreviousView);
			ShadowTemporalFilter.SetFloat("u_ClampBias", 0.001f);
			ShadowTemporalFilter.SetBool("u_DiffuseTemporal", false);
			ShadowTemporalFilter.SetBool("u_ShadowTemporal", true);

			ShadowTemporalFilter.SetFloat("u_MinimumMix", 0.0f);
			ShadowTemporalFilter.SetFloat("u_MaximumMix", ModifiedWorld ? 0.375250f : 0.95f);
			ShadowTemporalFilter.SetInteger("u_TemporalQuality", 1);
			ShadowTemporalFilter.SetBool("u_ReflectionTemporal", false);
			ShadowTemporalFilter.SetBool("u_ShouldFilterShadows", DenoiseSunShadows);

			ShadowTemporalFilter.SetMatrix4("u_VertInverseView", inv_view);
			ShadowTemporalFilter.SetMatrix4("u_VertInverseProjection", inv_projection);
			ShadowTemporalFilter.SetMatrix4("u_InverseView", inv_view);
			ShadowTemporalFilter.SetMatrix4("u_InverseProjection", inv_projection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ShadowRawTrace.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, PrevShadowTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(0));

			glActiveTexture(GL_TEXTURE8);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));
			
			glActiveTexture(GL_TEXTURE17);
			glBindTexture(GL_TEXTURE_2D, ShadowRawTrace.GetTexture(1));

			glActiveTexture(GL_TEXTURE22);
			glBindTexture(GL_TEXTURE_2D, PrevShadowTemporalFBO.GetTexture(1));
			
			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ShadowTemporalFBO.Unbind();
		}

		// Final shadow noise cleanup
		if (SoftShadows&&DenoiseSunShadows)
		{
			ShadowFiltered.Bind();
			ShadowFilter.Use();
			ShadowFilter.SetInteger("u_InputTexture", 0);
			ShadowFilter.SetInteger("u_PositionTexture", 1);
			ShadowFilter.SetInteger("u_NormalTexture", 2);
			ShadowFilter.SetInteger("u_IntersectionTransversals", 3);
			ShadowFilter.SetInteger("u_FrameCount", 4);
			ShadowFilter.SetVector2f("u_Dimensions", glm::vec2(ShadowFiltered.GetWidth(), ShadowFiltered.GetHeight()));
			ShadowFilter.SetMatrix4("u_VertInverseView", inv_view);
			ShadowFilter.SetMatrix4("u_VertInverseProjection", inv_projection);
			ShadowFilter.SetMatrix4("u_InverseView", inv_view);
			ShadowFilter.SetMatrix4("u_InverseProjection", inv_projection);
			ShadowFilter.SetFloat("u_ShadowFilterScale", ShadowFilterScale);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ShadowTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, ShadowRawTrace.GetTexture(1));

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D, ShadowTemporalFBO.GetTexture(1));

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}

		if (SSSSS && true && SSSSSStrength > 0.01f) {

			auto& FinalShadowfbo = SoftShadows ? (DenoiseSunShadows ? ShadowFiltered : ShadowTemporalFBO) : ShadowRawTrace;
			//SSSBlurPoisson
			
			// Poisson pass ->

			SSSBlurPoisson.Use();
			ShadowSSS.Bind();
			SSSBlurPoisson.SetInteger("u_Texture", 0);
			SSSBlurPoisson.SetInteger("u_BlockIDs", 1);
			SSSBlurPoisson.SetInteger("u_Depth", 2);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, FinalShadowfbo.GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			// Gaussian pass ->

			SSSBlur.Use();
			ShadowSSS2.Bind();
			SSSBlur.SetInteger("u_Texture", 0);
			SSSBlur.SetInteger("u_BlockIDs", 1);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ShadowSSS.GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}

		glUseProgram(0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// ---- REFLECTION TRACE ----


		{
			ReflectionTraceFBO.Bind();
			ReflectionTraceShader.Use();

			ReflectionTraceShader.SetInteger("u_PositionTexture", 0);
			ReflectionTraceShader.SetInteger("u_BlockIDTex", 1);
			ReflectionTraceShader.SetInteger("u_InitialTraceNormalTexture", 2);
			ReflectionTraceShader.SetInteger("u_BlockNormalTextures", 4);
			ReflectionTraceShader.SetInteger("u_BlockAlbedoTextures", 5);
			ReflectionTraceShader.SetInteger("u_BlockPBRTextures", 6);
			ReflectionTraceShader.SetInteger("u_Skymap", 7);
			ReflectionTraceShader.SetInteger("u_VoxelData", 8);
			ReflectionTraceShader.SetInteger("u_BlueNoiseTexture", 9);
			ReflectionTraceShader.SetInteger("u_DistanceFieldTexture", 10);
			ReflectionTraceShader.SetInteger("u_BlockEmissiveTextures", 11);
			ReflectionTraceShader.SetInteger("u_DiffuseSH", 14);
			ReflectionTraceShader.SetInteger("u_DiffuseCoCg", 15);
			ReflectionTraceShader.SetInteger("u_ShadowTrace", 16);

			// LPV ->
			ReflectionTraceShader.SetInteger("u_LPV", 20);
			ReflectionTraceShader.SetInteger("u_LPVBlocks", 21);
			ReflectionTraceShader.SetInteger("u_GBufferNormals", 23);
			ReflectionTraceShader.SetInteger("u_GBufferPBR", 24);

			ReflectionTraceShader.SetFloat("u_ReflectionTraceRes", ReflectionTraceResolution);
			ReflectionTraceShader.SetVector3f("u_SunDirection", SunDirection);
			ReflectionTraceShader.SetVector3f("u_MoonDirection", MoonDirection);
			ReflectionTraceShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
			ReflectionTraceShader.SetVector2f("u_Dimensions", glm::vec2(ReflectionTraceFBO.GetWidth(), ReflectionTraceFBO.GetHeight()));
			ReflectionTraceShader.SetVector2f("u_Halton", glm::vec2(GetTAAJitterSecondary(app.GetCurrentFrame())));
			ReflectionTraceShader.SetFloat("u_Time", glfwGetTime());
			ReflectionTraceShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());
			ReflectionTraceShader.SetBool("u_RoughReflections", RoughReflections);
			ReflectionTraceShader.SetBool("u_UseBlueNoise", USE_BLUE_NOISE_FOR_TRACING);
			ReflectionTraceShader.SetBool("u_ReflectPlayer", REFLECT_PLAYER);
			ReflectionTraceShader.SetBool("u_CloudReflections", CloudReflections);
			ReflectionTraceShader.SetBool("u_TemporalFilterReflections", TEMPORAL_SPEC);
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
			ReflectionTraceShader.SetInteger("u_PlayerSprite", 12);
			ReflectionTraceShader.SetInteger("u_IndirectAO", 17);
			ReflectionTraceShader.SetInteger("u_ProjectedClouds", 18);
			ReflectionTraceShader.SetInteger("u_ReflectionTraceLength", ReflectionTraceLength);
			ReflectionTraceShader.SetInteger("u_SPP", ReflectionSPP);
			ReflectionTraceShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
			ReflectionTraceShader.SetInteger("u_CurrentFrameMod128", app.GetCurrentFrame()%128);
			ReflectionTraceShader.SetFloat("u_SunStrengthModifier", SunStrengthModifier);
			ReflectionTraceShader.SetFloat("u_MoonStrengthModifier", MoonStrengthModifier);
			ReflectionTraceShader.SetMatrix4("u_VertInverseView", inv_view);
			ReflectionTraceShader.SetMatrix4("u_VertInverseProjection", inv_projection);
			ReflectionTraceShader.SetMatrix4("u_InverseView", inv_view);
			ReflectionTraceShader.SetMatrix4("u_InverseProjection", inv_projection);
			ReflectionTraceShader.SetMatrix4("u_View", MainCamera.GetViewMatrix());
			ReflectionTraceShader.SetMatrix4("u_Projection", MainCamera.GetProjectionMatrix());
			ReflectionTraceShader.SetBool("u_ReprojectToScreenSpace", ReprojectReflectionsToScreenSpace);
			ReflectionTraceShader.SetBool("CHECKERBOARD_SPEC_SPP", CHECKERBOARD_SPEC_SPP);
			ReflectionTraceShader.SetBool("TEMPORAL_SPEC", TEMPORAL_SPEC);
			ReflectionTraceShader.SetBool("u_LPVGI", ReflectionLPVGI);
			ReflectionTraceShader.SetBool("u_QualityLPVGI", ReflectionHighQualityLPVGI);
			ReflectionTraceShader.SetBool("u_RoughnessBias", ReflectionRoughnessBias&&DenoiseReflections);
			ReflectionTraceShader.SetBool("u_HandleLobeDeviation", ReflectionDenoiserDeviationHandling);
			ReflectionTraceShader.SetBool("u_DeriveFromDiffuseSH", DeriveReflectionsFromDiffuseSH);
			ReflectionTraceShader.SetBool("u_ScreenSpaceSkylightingValid", USE_SVGF);
			ReflectionTraceShader.SetBool("u_UseDecoupledGI", RefectionUseDecoupledGI);

			ReflectionTraceShader.SetInteger("u_LavaBlockID", BlockDatabase::GetBlockID("Lava"));


			//ReflectionTraceShader.BindUBOToBindingPoint("UBO_BlockData", 0);
			BlockDataStorageBuffer.Bind(0);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

			glActiveTexture(GL_TEXTURE2);

			if (RandomDebugVar && DOWNSAMPLE_GBUFFERS) {
				glBindTexture(GL_TEXTURE_2D, QuarterResGBuffer.GetTexture(1));
				
				if (app.GetCurrentFrame() % 30 == 0) {
					std::cout << "\nDEBUG VAR OUTPUT -> PASSED DOWNSAMPLED GBUFFER TO SPECULAR INDIRECT TRACE SHADER\n";
				}

			}

			else {
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));
			}

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE6);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

			glActiveTexture(GL_TEXTURE7);
			glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapSecondary.GetTexture());

			glActiveTexture(GL_TEXTURE8);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE9);
			glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE10);
			glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE11);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

			glActiveTexture(GL_TEXTURE12);
			glBindTexture(GL_TEXTURE_2D, PlayerSprite.GetTextureID());

			glActiveTexture(GL_TEXTURE14);
			glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture(0));

			glActiveTexture(GL_TEXTURE15);
			glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture(1));

			glActiveTexture(GL_TEXTURE16);
			glBindTexture(GL_TEXTURE_2D, SoftShadows ?(DenoiseSunShadows?ShadowFiltered.GetTexture():ShadowTemporalFBO.GetTexture()) : ShadowRawTrace.GetTexture());
			
			glActiveTexture(GL_TEXTURE17);
			glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture(3));

			glActiveTexture(GL_TEXTURE18);
			glBindTexture(GL_TEXTURE_2D, CloudProjection.GetTexture());

			glActiveTexture(GL_TEXTURE20);
			glBindTexture(GL_TEXTURE_3D, VoxelRT::Volumetrics::GetDensityVolume());

			glActiveTexture(GL_TEXTURE21);
			glBindTexture(GL_TEXTURE_3D, VoxelRT::Volumetrics::GetColorVolume());

			glActiveTexture(GL_TEXTURE23);
			glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(1));
			glActiveTexture(GL_TEXTURE24);
			glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(2));

			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, BlueNoise_SSBO.m_SSBO);
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, VoxelRT::Volumetrics::GetAverageColorSSBO());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionTraceFBO.Unbind();


			// denoise the hit distance 
			if (DENOISE_REFLECTION_HIT_DATA)
			{
				if (app.GetCurrentFrame() % 16 == 0) { std::cout << "\ndenoised hit data\n"; }


				// 1st spatial pass ->

				// non-conservative, step size = 2

				ReflectionHitDataDenoised.Bind();
				BilateralHitDist_1.Use();

				BilateralHitDist_1.SetInteger("u_PositionTexture", 0);
				BilateralHitDist_1.SetInteger("u_HitDist", 1);
				BilateralHitDist_1.SetInteger("u_Normals", 2);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture(1));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();

				ReflectionHitDataDenoised.Unbind();


				// 2nd spatial pass :-
				// step size = 4

				ReflectionTraceFBO.Bind();
				BilateralHitDist_2.Use();

				BilateralHitDist_2.SetInteger("u_PositionTexture", 0);
				BilateralHitDist_2.SetInteger("u_HitDist", 1);
				BilateralHitDist_2.SetInteger("u_Normals", 2);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, ReflectionHitDataDenoised.GetTexture());

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}
		}

		// Temporally filter it
		//if (RoughReflections)
		{
			ReflectionTemporalFBO.Bind();
			SpecularTemporalFilter.Use();

			SpecularTemporalFilter.SetInteger("u_CurrentColorTexture", 0);
			SpecularTemporalFilter.SetInteger("u_CurrentPositionTexture", 1);
			SpecularTemporalFilter.SetInteger("u_PreviousColorTexture", 2);
			SpecularTemporalFilter.SetInteger("u_PreviousFramePositionTexture", 3);
			SpecularTemporalFilter.SetInteger("u_PBRTex", 4);


			SpecularTemporalFilter.SetInteger("u_SpecularHitDist", 7);
			SpecularTemporalFilter.SetInteger("u_PrevSpecularHitDist", 8);
			SpecularTemporalFilter.SetInteger("u_NormalTexture", 9);
			SpecularTemporalFilter.SetInteger("u_PreviousNormalTexture", 10);
			SpecularTemporalFilter.SetInteger("u_EmissivityIntersectionMask", 11);
			SpecularTemporalFilter.SetInteger("u_TemporalHitDist", 12);


			SpecularTemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
			SpecularTemporalFilter.SetMatrix4("u_View", CurrentView);
			SpecularTemporalFilter.SetMatrix4("u_PrevProjection", PreviousProjection);
			SpecularTemporalFilter.SetMatrix4("u_PrevView", PreviousView);

			SpecularTemporalFilter.SetFloat("u_MinimumMix", 0.0f); // Brutal temporal filtering
			SpecularTemporalFilter.SetFloat("u_MaximumMix", 0.95f);
			SpecularTemporalFilter.SetInteger("u_TemporalQuality", 1);
			SpecularTemporalFilter.SetBool("u_ReflectionTemporal", true);
			SpecularTemporalFilter.SetBool("TEMPORAL_SPEC", TEMPORAL_SPEC);
			SpecularTemporalFilter.SetBool("u_FireflyRejection", ReflectionFireflyRejection);
			SpecularTemporalFilter.SetBool("u_AggressiveFireflyRejection", AggressiveFireflyRejection);
			SpecularTemporalFilter.SetBool("u_SmartClip", SmartReflectionClip);
			SpecularTemporalFilter.SetBool("u_RoughnessWeight", RoughReflections);
			SpecularTemporalFilter.SetBool("u_TemporallyStabializeHitDistance", TemporallyStabializeHitDistance);

			SpecularTemporalFilter.SetVector3f("u_PrevCameraPos", PreviousPosition);
			SpecularTemporalFilter.SetVector3f("u_CurrentCameraPos", MainCamera.GetPosition());

			SpecularTemporalFilter.SetMatrix4("u_VertInverseView", inv_view);
			SpecularTemporalFilter.SetMatrix4("u_VertInverseProjection", inv_projection);
			SpecularTemporalFilter.SetMatrix4("u_InverseView", inv_view);
			SpecularTemporalFilter.SetMatrix4("u_InverseProjection", inv_projection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, PrevReflectionTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(0));

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(2));



			glActiveTexture(GL_TEXTURE7);
			glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture(1));
			glActiveTexture(GL_TEXTURE8);
			glBindTexture(GL_TEXTURE_2D, PrevReflectionTraceFBO.GetTexture(1));

			glActiveTexture(GL_TEXTURE9);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE10);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(1));

			glActiveTexture(GL_TEXTURE11);
			glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture(2));

			glActiveTexture(GL_TEXTURE12);
			glBindTexture(GL_TEXTURE_2D, PrevReflectionTemporalFBO.GetTexture(2));


			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			ReflectionTemporalFBO.Unbind();
		}

		// Denoise reflection trace
		if (DenoiseReflections && RoughReflections)
		{
			// x pass
				
			{
				ReflectionDenoised_1.Bind();
				ReflectionDenoiser.Use();

				ReflectionDenoiser.SetInteger("u_InputTexture", 0);
				ReflectionDenoiser.SetInteger("u_PositionTexture", 1);
				ReflectionDenoiser.SetInteger("u_NormalTexture", 2);
				ReflectionDenoiser.SetInteger("u_PBRTex", 3);
				ReflectionDenoiser.SetInteger("u_BlockIDTex", 6);
				ReflectionDenoiser.SetInteger("u_BlockNormals", 11);
				ReflectionDenoiser.SetInteger("u_SpecularHitData", 10);
				ReflectionDenoiser.SetInteger("u_GBufferNormals", 13);
				ReflectionDenoiser.SetInteger("u_GBufferPBR", 14);
				ReflectionDenoiser.SetInteger("u_Frames", 15);
				ReflectionDenoiser.SetInteger("u_Step", 1);
				ReflectionDenoiser.SetBool("u_Dir", true);
				ReflectionDenoiser.SetBool("u_NormalMapAware", ReflectionNormalMapWeight);
				ReflectionDenoiser.SetBool("TEMPORAL_SPEC", TEMPORAL_SPEC);
				ReflectionDenoiser.SetBool("u_RoughnessBias", ReflectionRoughnessBias);
				ReflectionDenoiser.SetVector2f("u_Dimensions", glm::vec2(ReflectionDenoised_1.GetWidth(), ReflectionDenoised_1.GetHeight()));
				ReflectionDenoiser.SetMatrix4("u_VertInverseView", inv_view);
				ReflectionDenoiser.SetMatrix4("u_VertInverseProjection", inv_projection);
				ReflectionDenoiser.SetMatrix4("u_InverseView", inv_view);
				ReflectionDenoiser.SetMatrix4("u_InverseProjection", inv_projection);
				ReflectionDenoiser.SetMatrix4("u_PrevProjection", PreviousProjection);
				ReflectionDenoiser.SetMatrix4("u_PrevView", PreviousView);
				ReflectionDenoiser.SetMatrix4("u_View", MainCamera.GetViewMatrix());
				ReflectionDenoiser.SetMatrix4("u_Projection", MainCamera.GetProjectionMatrix());
				ReflectionDenoiser.SetBool("u_HandleLobeDeviation", ReflectionDenoiserDeviationHandling);
				ReflectionDenoiser.SetInteger("u_BlockIDTex", 6);
				ReflectionDenoiser.SetInteger("u_BlockPBRTexArray", 7);
				ReflectionDenoiser.SetInteger("u_ReflectionDenoisingRadiusBias", ReflectionDenoisingRadiusBias);
				ReflectionDenoiser.SetFloat("u_Time", glfwGetTime());
				ReflectionDenoiser.SetFloat("u_ResolutionScale", ReflectionSuperSampleResolution);
				ReflectionDenoiser.SetFloat("u_NormalMapWeightStrength", NormalMapWeightStrength);
				ReflectionDenoiser.SetBool("u_DeriveFromDiffuseSH", DeriveReflectionsFromDiffuseSH);
				ReflectionDenoiser.SetFloat("u_ReflectionDenoiserScale", ReflectionDenoiserScale);
				ReflectionDenoiser.SetBool("u_TemporalWeight", ReflectionTemporalWeight&&TEMPORAL_SPEC);
				ReflectionDenoiser.SetFloat("u_RoughnessNormalWeightBiasStrength", RoughnessNormalWeightBiasStrength);
				ReflectionDenoiser.SetBool("u_AmplifyReflectionTransversalWeight", AmplifyReflectionTransversalWeight);

				BlockDataStorageBuffer.Bind(0);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(2));


				glActiveTexture(GL_TEXTURE6);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

				glActiveTexture(GL_TEXTURE7);
				glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetPBRTextureArray());

				glActiveTexture(GL_TEXTURE10);

				if (TEMPORAL_SPEC&&TemporallyStabializeHitDistance) {
					glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture(2));
				}

				else {
					glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture(1));
				}

				glActiveTexture(GL_TEXTURE11);
				glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

				glActiveTexture(GL_TEXTURE13);
				glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(1));
				glActiveTexture(GL_TEXTURE14);
				glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(2));
				
				glActiveTexture(GL_TEXTURE15);
				glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture(1));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}

			// y pass 

			{
				ReflectionDenoised_2.Bind();
				ReflectionDenoiser.Use();

				ReflectionDenoiser.SetInteger("u_GBufferNormals", 13);
				ReflectionDenoiser.SetInteger("u_GBufferPBR", 14);
				ReflectionDenoiser.SetInteger("u_InputTexture", 0);
				ReflectionDenoiser.SetInteger("u_PositionTexture", 1);
				ReflectionDenoiser.SetInteger("u_NormalTexture", 2);
				ReflectionDenoiser.SetInteger("u_BlockIDTex", 6);
				ReflectionDenoiser.SetInteger("u_BlockPBRTexArray", 7);
				ReflectionDenoiser.SetInteger("u_SpecularHitData", 10);
				ReflectionDenoiser.SetInteger("u_BlockNormals", 11);
				ReflectionDenoiser.SetInteger("u_Frames", 15);
				ReflectionDenoiser.SetInteger("u_Step", 1);
				ReflectionDenoiser.SetBool("u_Dir", false);
				ReflectionDenoiser.SetBool("TEMPORAL_SPEC", TEMPORAL_SPEC);
				ReflectionDenoiser.SetVector2f("u_Dimensions", glm::vec2(ReflectionDenoised_1.GetWidth(), ReflectionDenoised_1.GetHeight()));
				ReflectionDenoiser.SetMatrix4("u_VertInverseView", inv_view);
				ReflectionDenoiser.SetMatrix4("u_VertInverseProjection", inv_projection);
				ReflectionDenoiser.SetMatrix4("u_InverseView", inv_view);
				ReflectionDenoiser.SetMatrix4("u_InverseProjection", inv_projection);
				ReflectionDenoiser.SetMatrix4("u_PrevProjection", PreviousProjection);
				ReflectionDenoiser.SetMatrix4("u_PrevView", PreviousView);
				ReflectionDenoiser.SetFloat("u_Time", glfwGetTime());
				ReflectionDenoiser.SetFloat("u_ResolutionScale", ReflectionSuperSampleResolution);
				ReflectionDenoiser.SetBool("u_NormalMapAware", ReflectionNormalMapWeight);
				ReflectionDenoiser.SetBool("u_TemporalWeight", ReflectionTemporalWeight);
				ReflectionDenoiser.SetInteger("u_ReflectionDenoisingRadiusBias", ReflectionDenoisingRadiusBias);
				ReflectionDenoiser.SetBool("u_HandleLobeDeviation", ReflectionDenoiserDeviationHandling);
				ReflectionDenoiser.SetBool("u_DeriveFromDiffuseSH", DeriveReflectionsFromDiffuseSH);
				ReflectionDenoiser.SetFloat("u_NormalMapWeightStrength", NormalMapWeightStrength);
				ReflectionDenoiser.SetFloat("u_RoughnessNormalWeightBiasStrength", RoughnessNormalWeightBiasStrength);
				ReflectionDenoiser.SetFloat("u_ReflectionDenoiserScale", ReflectionDenoiserScale);
				ReflectionDenoiser.SetBool("u_TemporalWeight", ReflectionTemporalWeight&& TEMPORAL_SPEC);
				ReflectionDenoiser.SetBool("u_AmplifyReflectionTransversalWeight", AmplifyReflectionTransversalWeight);

				BlockDataStorageBuffer.Bind(0);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, ReflectionDenoised_1.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));


				glActiveTexture(GL_TEXTURE6);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

				glActiveTexture(GL_TEXTURE7);
				glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetPBRTextureArray());

				glActiveTexture(GL_TEXTURE10);

				if (TEMPORAL_SPEC&&TemporallyStabializeHitDistance) {
					glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture(2));
				}

				else {
					glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture(1));
				}

				glActiveTexture(GL_TEXTURE11);
				glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

				glActiveTexture(GL_TEXTURE13);
				glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(1));
				glActiveTexture(GL_TEXTURE14);
				glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(2));

				glActiveTexture(GL_TEXTURE15);
				glBindTexture(GL_TEXTURE_2D, ReflectionTemporalFBO.GetTexture(1));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}
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
			RTAOShader.SetInteger("u_BlockIDTexture", 5);
			RTAOShader.SetFloat("u_Time", glfwGetTime());
			RTAOShader.SetMatrix4("u_VertInverseView", inv_view);
			RTAOShader.SetMatrix4("u_VertInverseProjection", inv_projection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

			BlockDataStorageBuffer.Bind(0);

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
			MainTemporalFilter.SetBool("u_DiffuseTemporal", false);
			MainTemporalFilter.SetBool("u_ShadowTemporal", false);

			MainTemporalFilter.SetFloat("u_MinimumMix", 0.25f);
			MainTemporalFilter.SetFloat("u_MaximumMix", 0.975f);
			MainTemporalFilter.SetInteger("u_TemporalQuality", 0);
			MainTemporalFilter.SetBool("u_ReflectionTemporal", false);

			MainTemporalFilter.SetMatrix4("u_VertInverseView", inv_view);
			MainTemporalFilter.SetMatrix4("u_VertInverseProjection", inv_projection);
			MainTemporalFilter.SetMatrix4("u_InverseView", inv_view);
			MainTemporalFilter.SetMatrix4("u_InverseProjection", inv_projection);
			MainTemporalFilter.SetFloat("u_ClampBias", 0.02f);
			MainTemporalFilter.SetInteger("u_NormalTexture", 8);

			glActiveTexture(GL_TEXTURE8);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, RTAO_FBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, PrevRTAOTemporalFBO.GetTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(0));

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			RTAOTemporalFBO.Unbind();
		}

		// ---- RENDER CLOUDS ----
		GLuint CloudData = 0;
		if (CloudsEnabled)
		{
			bool UpdateCloudProjection = (CloudReflections) || app.GetCurrentFrame() % 2 == 0;

			glm::vec3 CelestialDirections[2] = { SunDirection, MoonDirection };

			CloudData = Clouds::CloudRenderer::Update(MainCamera.GetProjectionMatrix(), MainCamera.GetViewMatrix(), PreviousProjection,
				PreviousView, CurrentPosition,
				PreviousPosition, VAO, StrongerLightDirection, BluenoiseTexture.GetTextureID(),
				PADDED_WIDTH, PADDED_HEIGHT, app.GetCurrentFrame(), SkymapSecondary.GetTexture(), InitialTraceFBO->GetTexture(0), PreviousPosition, InitialTraceFBOPrev->GetTexture(0),
				CloudModifiers, ClampCloudTemporal, glm::vec3(CloudDetailScale,CloudDetailWeightEnabled?1.0f:0.0f,CloudErosionWeightExponent), 
				CloudTimeScale, CurlNoiseOffset, CirrusStrength,CirrusScale, CloudStepCount, CloudCheckerStepCount, sun_visibility, CloudDetailFBMPower, 
				CloudLODLighting, CloudForceSupersample, CloudForceSupersampleRes, CloudSpatialUpscale, CloudAmbientDensityMultiplier, CloudProjection.GetTexture(0), UpdateCloudProjection, CloudThiccness, CloudDetailScale, false, CelestialDirections, SkymapMain.GetTexture(), glm::vec2(SunStrengthModifier, MoonStrengthModifier));
			


			Clouds::CloudRenderer::SetCoverage(CloudCoverage);
			Clouds::CloudRenderer::SetResolution(CloudResolution);
		}

		// ---- COLOR PASS ----

		BlockDatabase::SetTextureArraysFilteringNearest();

		ReflectionProjection = PreviousProjection;
		ReflectionView = PreviousView;

		ColoredFBO.Bind();

		ColorShader.Use();
		ColorShader.SetInteger("u_DiffuseTexture", 0);
		ColorShader.SetInteger("u_NormalTexture", 1);
		ColorShader.SetInteger("u_InitialTracePositionTexture", 2);
		ColorShader.SetInteger("u_BlockIDTexture", 3);
		ColorShader.SetInteger("u_BlockAlbedoTextures", 4);
		ColorShader.SetInteger("u_BlockNormalTextures", 5);
		ColorShader.SetInteger("u_BlockPBRTextures", 6);
		ColorShader.SetInteger("u_Skybox", 7);
		ColorShader.SetInteger("u_ShadowTexture", 8);
		ColorShader.SetInteger("u_BlueNoiseTextures", 9);
		ColorShader.SetInteger("u_BlockEmissiveTextures", 11);
		ColorShader.SetInteger("u_CloudData", 12);

		// indirect diffuse sh
		ColorShader.SetInteger("u_DiffuseSHData1", 14);
		ColorShader.SetInteger("u_DiffuseSHData2", 15);

		ColorShader.SetInteger("u_DiffuseSHy", 14);
		ColorShader.SetInteger("u_DiffuseCoCg", 15);

		ColorShader.SetInteger("u_ReflectionSHData", 16);
		ColorShader.SetInteger("u_SpecularIndirect", 16);
		ColorShader.SetInteger("u_ReflectionCoCgData", 17);
		ColorShader.SetInteger("u_VXAO", 18);
		ColorShader.SetInteger("u_HighResBL", 19);
		
		ColorShader.SetInteger("u_LPVLightLevel", 20);
		ColorShader.SetInteger("u_LPVColorData", 21);
		ColorShader.SetInteger("u_BlueNoiseHighRes", 23);
		ColorShader.SetInteger("u_LPVDebugState", LPVDebugState);

		ColorShader.SetInteger("u_ContactHardeningShadows", SoftShadows);
		ColorShader.SetMatrix4("u_InverseView", inv_view);
		ColorShader.SetMatrix4("u_InverseProjection", inv_projection);
		ColorShader.SetMatrix4("u_View", MainCamera.GetViewMatrix());
		ColorShader.SetMatrix4("u_Projection", MainCamera.GetProjectionMatrix());
		ColorShader.SetMatrix4("u_ShadowProjection", ShadowProjection);
		ColorShader.SetMatrix4("u_ShadowView", ShadowView);
		ColorShader.SetMatrix4("u_ReflectionProjection", ReflectionProjection);
		ColorShader.SetMatrix4("u_ReflectionView", ReflectionView);
		ColorShader.SetVector2f("u_InitialTraceResolution", glm::vec2(floor(PADDED_WIDTH * InitialTraceResolution), floor(PADDED_HEIGHT * InitialTraceResolution)));
		ColorShader.SetVector3f("u_SunDirection", glm::normalize(SunDirection));
		ColorShader.SetVector3f("u_MoonDirection", glm::normalize(MoonDirection));
		ColorShader.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
		ColorShader.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());

		ColorShader.SetVector3f("u_StrongerLightDirectionVert", StrongerLightDirection);
		ColorShader.SetVector3f("u_SunDirectionVert", SunDirection);
		ColorShader.SetVector3f("u_MoonDirectionVert", MoonDirection);

		ColorShader.SetFloat("u_Time", glfwGetTime());
		ColorShader.SetFloat("u_GrassblockAlbedoID", BlockDatabase::GetBlockTexture("Grass", BlockDatabase::BlockFaceType::Front));
		ColorShader.SetFloat("u_SunStrengthModifier", SunStrengthModifier);
		ColorShader.SetFloat("u_MoonStrengthModifier", MoonStrengthModifier);
		ColorShader.SetFloat("u_TextureDesatAmount", TextureDesatAmount);
		ColorShader.SetFloat("u_SubsurfaceScatterStrength", SSSSSStrength);
		ColorShader.SetBool("u_CloudsEnabled", CloudsEnabled);
		ColorShader.SetBool("u_VXAOCutoff", VXAO_CUTOFF);
		ColorShader.SetBool("u_POM", POM);
		ColorShader.SetFloat("u_POMHeight", POMHeight);
		ColorShader.SetBool("u_HighQualityPOM", HighQualityPOM);
		ColorShader.SetBool("u_DitherPOM", DitherPOM);

		ColorShader.SetBool("u_RemoveTiling", RemoveTiling);
		ColorShader.SetBool("u_RTAO", RTAO);
		ColorShader.SetBool("u_AmplifyNormalMap", AmplifyNormalMap);
		ColorShader.SetBool("u_DitherPOM", DitherPOM);
		ColorShader.SetBool("u_DoVXAO", VXAO);
		ColorShader.SetBool("u_SVGFEnabled", USE_SVGF);
		ColorShader.SetBool("u_DEBUGDiffuseGI", DEBUGTraceLevel==1);
		ColorShader.SetBool("u_DEBUGSpecGI", DEBUGTraceLevel==2);
		ColorShader.SetBool("u_DEBUGShadows", DEBUGTraceLevel==3);
		ColorShader.SetBool("u_ShouldDitherUpscale", DITHER_SPATIAL_UPSCALE);
		ColorShader.SetBool("u_UseDFG", UseEnvironmentBRDF);
		ColorShader.SetBool("u_SSSSS", SSSSS && SSSSSStrength > 0.01f);
		ColorShader.SetBool("u_NebulaCelestialColor", NebulaCelestialColor);
		ColorShader.SetBool("u_InferSpecularDetailSpatially", InferSpecularDetailSpatially);
		ColorShader.SetBool("u_CloudCatmullRomUpsampling", CloudFinalCatmullromUpsample);
		ColorShader.SetVector2f("u_Dimensions", glm::vec2(PADDED_WIDTH, PADDED_HEIGHT));
		ColorShader.SetVector2f("u_Halton", GetTAAJitterSecondary(app.GetCurrentFrame()));
		ColorShader.SetMatrix4("u_InverseView", inv_view);
		ColorShader.SetMatrix4("u_InverseProjection", inv_projection);
		ColorShader.SetInteger("u_DebugTexture", 22);
		ColorShader.SetInteger("u_CloudProjectionVert", 28);
		ColorShader.SetInteger("u_CloudProjection", 28);
		ColorShader.SetInteger("u_SSSShadowMap", 24);
		ColorShader.SetVector4f("u_FocusedOnBlock", HighlightFocusedBlock ? glm::vec4(FocusedOnBlock) : glm::vec4(-1.0f));

		ColorShader.SetInteger("u_LavaBlockID", BlockDatabase::GetBlockID("Lava"));
		ColorShader.SetInteger("u_LavaTextures[0]", 26);
		ColorShader.SetInteger("u_LavaTextures[1]", 27);


		ColorShader.SetInteger("u_NebulaLowRes", 25);



		ColorShader.SetInteger("u_GBufferAlbedos", 29);
		ColorShader.SetInteger("u_GBufferNormals", 30);
		ColorShader.SetInteger("u_GBufferPBR", 31);
		ColorShader.SetInteger("u_GBufferAO", 10);


		ColorShader.SetFloat("u_NebulaStrength", NebulaStrength);



		

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapMain.GetTexture());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, SoftShadows ? (DenoiseSunShadows ? ShadowFiltered.GetTexture() : ShadowTemporalFBO.GetTexture()) : ShadowRawTrace.GetTexture());

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D_ARRAY, BlueNoise.GetTextureArray());

		//glActiveTexture(GL_TEXTURE10);
		//glBindTexture(GL_TEXTURE_2D, RoughReflections ? (DenoiseReflections ? ReflectionDenoised_2.GetTexture() : ReflectionTemporalFBO.GetTexture()) : ReflectionCheckerReconstructed.GetTexture());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());

		glActiveTexture(GL_TEXTURE12);
		glBindTexture(GL_TEXTURE_2D, CloudData);

		glActiveTexture(GL_TEXTURE14);
		glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture(0));
		//glBindTexture(GL_TEXTURE_2D, DiffuseTraceFBO.GetTexture(0));

		glActiveTexture(GL_TEXTURE15);
		glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture(1));
		//glBindTexture(GL_TEXTURE_2D, DiffuseTraceFBO.GetTexture(1));

		glActiveTexture(GL_TEXTURE16);
		glBindTexture(GL_TEXTURE_2D, RoughReflections ? (DenoiseReflections ? ReflectionDenoised_2.GetTexture() : ReflectionTemporalFBO.GetTexture()) : ReflectionTemporalFBO.GetTexture());

		glActiveTexture(GL_TEXTURE18);
		glBindTexture(GL_TEXTURE_2D, DiffuseDenoiseFBO.GetTexture(3));

		glActiveTexture(GL_TEXTURE19);
		glBindTexture(GL_TEXTURE_2D, BluenoiseHighResTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE20);
		glBindTexture(GL_TEXTURE_3D, VoxelRT::Volumetrics::GetDensityVolume());

		glActiveTexture(GL_TEXTURE21);
		glBindTexture(GL_TEXTURE_3D, VoxelRT::Volumetrics::GetColorVolume());
		

	    
		glActiveTexture(GL_TEXTURE22);
		glBindTexture(GL_TEXTURE_2D, CloudProjection.GetTexture());
		//glBindTexture(GL_TEXTURE_2D, ReflectionTraceFBO.GetTexture(2));
		//glBindTexture(GL_TEXTURE_2D, DiffuseRawTraceFBO.GetTexture(4));

		glActiveTexture(GL_TEXTURE23);
		glBindTexture(GL_TEXTURE_2D, BluenoiseHighResTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE24);
		glBindTexture(GL_TEXTURE_2D, ShadowSSS2.GetTexture());

		glActiveTexture(GL_TEXTURE25);
		glBindTexture(GL_TEXTURE_CUBE_MAP, NightSkyMapLowRes.GetID());


		glActiveTexture(GL_TEXTURE28);
		glBindTexture(GL_TEXTURE_2D, CloudProjection.GetTexture());



		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE29);
		glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(0));

		glActiveTexture(GL_TEXTURE30);
		glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(1));

		glActiveTexture(GL_TEXTURE31);
		glBindTexture(GL_TEXTURE_2D, GeneratedGBuffer.GetTexture(2));



		BlockDataStorageBuffer.Bind(0);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, VoxelRT::Volumetrics::GetAverageColorSSBO());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		ColoredFBO.Unbind();

		BlockDatabase::SetTextureArraysFilteringNearest();



		// ----- TAA ----- //



		TAAFBO.Bind();

		TemporalAAShader.Use();
		TemporalAAShader.SetInteger("u_CurrentColorTexture", 0);
		TemporalAAShader.SetInteger("u_PositionTexture", 1);
		TemporalAAShader.SetInteger("u_PreviousColorTexture", 2);
		TemporalAAShader.SetInteger("u_PreviousPositionTexture", 3);
		TemporalAAShader.SetBool("u_Enabled", TAA);
		TemporalAAShader.SetBool("u_BlockModified", ModifiedWorld);
		TemporalAAShader.SetBool("u_DepthWeight", TAADepthWeight);
		TemporalAAShader.SetFloat("u_DepthExponentMultiplier", TAADepthWeightExp);

		TemporalAAShader.SetMatrix4("u_PrevProjection", PreviousProjection);
		TemporalAAShader.SetMatrix4("u_PrevView", PreviousView);
		TemporalAAShader.SetMatrix4("u_InversePrevProjection", glm::inverse(PreviousProjection));
		TemporalAAShader.SetMatrix4("u_InversePrevView", glm::inverse(PreviousView));

		TemporalAAShader.SetMatrix4("u_VertInverseView", inv_view);
		TemporalAAShader.SetMatrix4("u_VertInverseProjection", inv_projection);
		TemporalAAShader.SetMatrix4("u_InverseView", inv_view);
		TemporalAAShader.SetMatrix4("u_InverseProjection", inv_projection);
		TemporalAAShader.SetMatrix4("u_View", MainCamera.GetViewMatrix());
		TemporalAAShader.SetMatrix4("u_Projection", MainCamera.GetProjectionMatrix());

		TemporalAAShader.SetVector3f("u_CameraHistory[0]", MainCamera.GetPosition());
		TemporalAAShader.SetVector3f("u_CameraHistory[1]", PreviousPosition);
		TemporalAAShader.SetVector2f("u_Halton", GetTAAJitter(app.GetCurrentFrame(), glm::vec2(TAAFBO.GetWidth(), TAAFBO.GetHeight())));

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetColorTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevTAAFBO.GetTexture());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBOPrev->GetTexture(3));

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
			SSAOShader.SetFloat("u_SSAOStrength", SSAOStrength);
			SSAOShader.SetMatrix4("u_VertInverseView", inv_view);
			SSAOShader.SetMatrix4("u_VertInverseProjection", inv_projection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

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

		GLuint BrightTex = 0;

		if (Bloom)
		{
			BloomRenderer::RenderBloom(BloomFBO, BloomFBOAlternate, GeneratedGBuffer.GetTexture(0), ColoredFBO.GetEmissiveMask(), false, BrightTex, BloomWide);
		}



		// ---- Auto Exposure ----



		// WIP  !

		float ComputedExposure = 3.25f;
		
		if (AutoExposure)
		{
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, AutoExposurePDF);
			glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(GLuint) * 256, &Zero256UINT);
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

			AutoExposureComputePDF.Use();
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, AutoExposurePDF);
			AutoExposureComputePDF.SetInteger("u_ColorTexture", 0);
			AutoExposureComputePDF.SetVector2f("u_Resolution", glm::vec2(ColoredFBO.GetWidth(), ColoredFBO.GetHeight()));
			AutoExposureComputePDF.SetVector2f("u_InverseResolution", glm::vec2(1.0f/ColoredFBO.GetWidth(), 1.0f/ColoredFBO.GetHeight()));
			AutoExposureComputePDF.SetFloat("u_DeltaTime", DeltaTime);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetColorTexture());

			glDispatchCompute(2, 2, 1);
			glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

			glUseProgram(0);

			ComputeAutoExposure.Use();
			ComputeAutoExposure.SetVector2f("u_Resolution", glm::vec2(ColoredFBO.GetWidth(), ColoredFBO.GetHeight()));
			ComputeAutoExposure.SetVector2f("u_InverseResolution", glm::vec2(1.0f / ColoredFBO.GetWidth(), 1.0f / ColoredFBO.GetHeight()));
			ComputeAutoExposure.SetFloat("u_DeltaTime", DeltaTime);

			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, AutoExposureSSBO);
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, AutoExposurePDF);

			glDispatchCompute(1, 1, 1);
			glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

			glUseProgram(0);

			float DownloadedExposure = 0.0f;
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, AutoExposureSSBO);
			glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(float), &DownloadedExposure);;
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
			

			// Compute final exposure ->
			float AE_BoundsMin = 3.0f;
			float AE_BoundsMax = 20.0f;
			float ExposureMultiplier = 5.0f;

			glm::vec3 AE_NormalizedSunDir = glm::normalize(SunDirection);
			float AE_SunVisibility = glm::clamp(glm::dot(AE_NormalizedSunDir, glm::vec3(0.0f, 1.0f, 0.0f)) + 0.05f, 0.0f, 0.1f) * 12.0;
			AE_SunVisibility = 1.0f - AE_SunVisibility;

			AE_BoundsMin = glm::mix(AE_BoundsMin, 2.35f, AE_SunVisibility);
			AE_BoundsMax = glm::mix(AE_BoundsMax, 3.01f, AE_SunVisibility);
			ExposureMultiplier = glm::mix(ExposureMultiplier, 1.0f, AE_SunVisibility);
			
			ComputedExposure = glm::clamp((DownloadedExposure * ExposureMultiplier), AE_BoundsMin, AE_BoundsMax);
			
			if (app.GetCurrentFrame() % 6 == 0) {
				std::cout << "\n AUTO EXPOSURE COMPUTED : " << ComputedExposure << "\n";
			}

			
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
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

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
			BilateralBlur.SetMatrix4("u_VertInverseView", inv_view);
			BilateralBlur.SetMatrix4("u_VertInverseProjection", inv_projection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, VolumetricFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

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
			BilateralBlur.SetMatrix4("u_VertInverseView", inv_view);
			BilateralBlur.SetMatrix4("u_VertInverseProjection", inv_projection);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, BlurredVolumetricFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			VolumetricFBO.Unbind();
		}

		// ACTUAL WORLD SPACE VOLUMETRICS //

		if (PointVolumetricsToggled && PointVolumetricStrength > 0.01f)
		{
			PointVolumetrics.Use();
			VolumetricsCompute.Bind();

			PointVolumetrics.SetMatrix4("u_ProjectionMatrix", MainCamera.GetProjectionMatrix());
			PointVolumetrics.SetMatrix4("u_ViewMatrix", MainCamera.GetViewMatrix());
			PointVolumetrics.SetMatrix4("u_InverseView", inv_view);
			PointVolumetrics.SetMatrix4("u_InverseProjection", inv_projection);
			PointVolumetrics.SetVector3f("u_ViewerPosition", MainCamera.GetPosition());

			PointVolumetrics.SetInteger("u_ParticipatingMedia", 0);
			PointVolumetrics.SetInteger("u_BlueNoise", 1);
			PointVolumetrics.SetInteger("u_LinearDepthTexture", 2);
			PointVolumetrics.SetInteger("u_VolumetricDensityData", 3);
			PointVolumetrics.SetInteger("u_VolumetricColorDataSampler", 4);
			PointVolumetrics.SetInteger("u_LightCount", world->LightChunkData.size());

			PointVolumetrics.SetFloat("u_Time", glfwGetTime());
			PointVolumetrics.SetFloat("u_Strength", PointVolumetricStrength);

			PointVolumetrics.SetBool("u_Colored", ColoredPointVolumetrics);
			PointVolumetrics.SetBool("u_UseBayer", PointVolumetricsBayer);
			PointVolumetrics.SetBool("u_UsePerlinNoiseForOD", PointVolPerlinOD);
			PointVolumetrics.SetBool("u_FractalSimplexOD", false); // FALSE
			PointVolumetrics.SetBool("u_PointVolTriquadraticDensityInterpolation", PointVolTriquadraticDensityInterpolation); 
			PointVolumetrics.SetBool("u_GroundTruthColorInterpolation", PointVolGroundTruthColorInterpolation); 

			PointVolumetrics.SetVector2f("u_Dimensions", glm::vec2(VolumetricsCompute.GetWidth(), VolumetricsCompute.GetHeight()));

			//glBindImageTexture(1, VoxelRT::Volumetrics::GetColorVolume(), 0, true, 0, GL_READ_ONLY, GL_R8UI);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_3D, VoxelRT::Volumetrics::GetDensityVolume());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_3D, VoxelRT::Volumetrics::GetColorVolume());

			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, VoxelRT::Volumetrics::GetAverageColorSSBO());
			//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, world->m_LightPositionSSBO);

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			// Denoise 

			if (DenoisePointVol) {
				VolumetricsComputeBlurred.Bind();
				VolumetricsDenoiser.Use();

				VolumetricsDenoiser.SetInteger("u_InputPointVolumetrics", 0);
				VolumetricsDenoiser.SetInteger("u_DepthTexture", 1);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, VolumetricsCompute.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}
			//// Blur 2
			//
			//Gaussian5TapOptimized.Use();
			//VolumetricsCompute.Bind();
			//
			//Gaussian5TapOptimized.SetInteger("u_Texture", 0);
			//
			//glActiveTexture(GL_TEXTURE0);
			//glBindTexture(GL_TEXTURE_2D, VolumetricsComputeBlurred.GetTexture());
			//
			//VAO.Bind();
			//glDrawArrays(GL_TRIANGLES, 0, 6);
			//VAO.Unbind();
		}

		

		// Combine bloom 


		if (Bloom) {
			BloomCombined.Bind();
			CombineBloom.Use();
			CombineBloom.SetInteger("u_BloomMips[0]", 0);
			CombineBloom.SetInteger("u_BloomMips[1]", 1);
			CombineBloom.SetInteger("u_BloomMips[2]", 2);
			CombineBloom.SetInteger("u_BloomMips[3]", 3);
			CombineBloom.SetInteger("u_BloomMips[4]", 4);
			CombineBloom.SetInteger("u_BloomBrightTexture", 5);
			CombineBloom.SetBool("u_HQBloomUpscale", HQBloomUpscale);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[0]);

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[1]);

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[2]);

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[3]);

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[4]);

			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_2D, BrightTex);

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}


		if (DoDiffractionSpikes) {
			DiffractionSpikesShader.Use();
			DiffractionSpikes.Bind();

			DiffractionSpikesShader.SetInteger("u_Input", 0);
			DiffractionSpikesShader.SetInteger("u_KernelSize", DiffractionKernel);
			DiffractionSpikesShader.SetFloat("u_Scale", DiffractionScaler);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[0]);
			//glBindTexture(GL_TEXTURE_2D, BloomCombined.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			if (BlurDiffraction) {

				Gaussian5TapOptimized.Use();
				DiffractionSpikesDenoised.Bind();

				Gaussian5TapOptimized.SetInteger("u_Texture", 0);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, DiffractionSpikes.GetTexture());

				VAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				VAO.Unbind();
			}
			
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
		PostProcessingShader.SetInteger("u_VolumetricTexture", 10);
		PostProcessingShader.SetInteger("u_RTAOTexture", 11);
		PostProcessingShader.SetInteger("u_NormalTexture", 12);
		PostProcessingShader.SetInteger("u_CloudData", 15);
		PostProcessingShader.SetInteger("u_VolumetricsCompute", 16);
		PostProcessingShader.SetInteger("u_BlockIDTexture", 17);
		PostProcessingShader.SetInteger("u_BloomCombined", 5);
		PostProcessingShader.SetInteger("u_DiffractionSpikes", 6);
		PostProcessingShader.SetInteger("u_BloomMip", 7);
		PostProcessingShader.SetInteger("u_DiffractionSpikesEnabled", DoDiffractionSpikes);
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
		PostProcessingShader.SetBool("u_HejlBurgess", HejlBurgessTonemap);
		PostProcessingShader.SetBool("u_Fucktard", Fucktard);
		//PostProcessingShader.SetBool("u_PurkinjeEffect", PurkinjeEffect);
		PostProcessingShader.SetBool("u_ExponentialFog", ExponentialFog);
		PostProcessingShader.SetBool("u_AutoExposure", AutoExposure);
		PostProcessingShader.SetBool("u_PointVolumetricsToggled", PointVolumetricsToggled&& PointVolumetricStrength > 0.01f);
		PostProcessingShader.SetFloat("u_LensFlareIntensity", LensFlareIntensity);
		PostProcessingShader.SetFloat("u_Exposure", ComputedExposure * ExposureMultiplier);
		PostProcessingShader.SetFloat("u_ChromaticAberrationStrength", ChromaticAberrationStrength);
		PostProcessingShader.SetFloat("u_GodRaysStrength", GodRaysStrength);
		PostProcessingShader.SetFloat("u_PurkingeEffectStrength", PurkingeEffectStrength);
		PostProcessingShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());
		PostProcessingShader.SetInteger("u_CurrentFrame", app.GetCurrentFrame());

		glm::mat4 CuberotationMap = glm::rotate(glm::mat4(1.0f), (float)glm::radians(glfwGetTime() * 2.0f), glm::vec3(0.3f, 0.1f, 1.0f));

		PostProcessingShader.SetInteger("u_NightSkymap", 21);
		PostProcessingShader.SetMatrix4("u_RotationMatrix", CuberotationMap);



		bool lff = LensFlare && StrongerLightDirection == SunDirection && LensFlareIntensity > 0.001f;



		PostProcessingShader.SetBool("u_LensFlare", lff);


		// Set the bloom mips
		
		PostProcessingShader.SetInteger("u_LensDirtOverlay", 24);
		PostProcessingShader.SetInteger("u_SmallestBloomMip", 25);
		PostProcessingShader.SetInteger("u_ShadowTexture", 9);
		PostProcessingShader.SetInteger("u_Clouds", 20);

		PostProcessingShader.SetMatrix4("u_VertInverseView", inv_view);
		PostProcessingShader.SetMatrix4("u_VertInverseProjection", inv_projection);
		PostProcessingShader.SetMatrix4("u_InverseView", inv_view);
		PostProcessingShader.SetMatrix4("u_InverseProjection", inv_projection);

		PostProcessingShader.SetVector3f("u_VertSunDir", SunDirection);
		PostProcessingShader.SetBool("u_ComputePlayerShadow", lff);
		PostProcessingShader.SetBool("u_FilmGrain", FilmGrainStrength>0.001f);

		PostProcessingShader.SetInteger("u_DistanceFieldTexture", 18);
		PostProcessingShader.SetInteger("u_VoxelVolume", 19);
		PostProcessingShader.SetInteger("u_Sky", 26);

		PostProcessingShader.SetFloat("u_Time", glfwGetTime());
		PostProcessingShader.SetFloat("u_FilmGrainStrength", FilmGrainStrength);
		PostProcessingShader.SetFloat("u_ExposureMultiplier", ExposureMultiplier);
		PostProcessingShader.SetFloat("u_NebulaStrength", NebulaStrength);
		PostProcessingShader.SetFloat("u_BloomStrength", BloomStrength);
		PostProcessingShader.SetFloat("u_DiffractionStrength", DiffractionStrength);

		PostProcessingShader.SetBool("u_LensDirt", LensDirt);
		PostProcessingShader.SetFloat("u_LensDirtStrength", LensDirtStrength);
		PostProcessingShader.SetBool("u_HQBloomUpscale", HQBloomUpscale);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TAAFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, BluenoiseTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, SSAOBlurred.GetTexture());
		
		
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, BloomCombined.GetTexture());
		
		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, BlurDiffraction ? DiffractionSpikesDenoised.GetTexture() : DiffractionSpikes.GetTexture());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[0]);


		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, SoftShadows ? (DenoiseSunShadows ? ShadowFiltered.GetTexture() : ShadowTemporalFBO.GetTexture()) : ShadowRawTrace.GetTexture());

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, VolumetricFBO.GetTexture());

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D, RTAOTemporalFBO.GetTexture());

		glActiveTexture(GL_TEXTURE12);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

		glActiveTexture(GL_TEXTURE13);
		glBindTexture(GL_TEXTURE_2D, ColoredFBO.GetEmissiveMask());

		glActiveTexture(GL_TEXTURE15);
		glBindTexture(GL_TEXTURE_2D, CloudData);

		glActiveTexture(GL_TEXTURE16);
		glBindTexture(GL_TEXTURE_2D, DenoisePointVol ?  VolumetricsComputeBlurred.GetTexture() : VolumetricsCompute.GetTexture());

		glActiveTexture(GL_TEXTURE17);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

		glActiveTexture(GL_TEXTURE18);
		glBindTexture(GL_TEXTURE_3D, world->m_DistanceFieldTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE19);
		glBindTexture(GL_TEXTURE_3D, world->m_DataTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE20);
		glBindTexture(GL_TEXTURE_2D, CloudData);

		glActiveTexture(GL_TEXTURE21);
		glBindTexture(GL_TEXTURE_CUBE_MAP, NightSkyMap.GetID());
		
			
		glActiveTexture(GL_TEXTURE24);
		glBindTexture(GL_TEXTURE_2D, LensDirtTexture.GetTextureID());

		glActiveTexture(GL_TEXTURE25);
		glBindTexture(GL_TEXTURE_2D, BloomFBO.m_Mips[4]);
		
		glActiveTexture(GL_TEXTURE26);
		glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapMain.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		PostProcessingFBO.Unbind();

		// DOF ->


		if (DOF) {

			// Downsample ->
			BicubicDownsampler.Use();
			DOFFBOInput.Bind();
			
			BicubicDownsampler.SetInteger("u_Texture", 0);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, PostProcessingFBO.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();


			// Blur ->
			DOFShader.Use();
			DOFFBO.Bind();

			DOFShader.SetInteger("u_DepthTexture", 0);
			DOFShader.SetInteger("u_InputTexture", 1);
			DOFShader.SetVector2f("u_Dimensions", glm::vec2(DOFFBO.GetWidth(), DOFFBO.GetHeight()));
			DOFShader.SetFloat("u_TexelScale", MainCamera.GetProjectionMatrix()[1][1] / 1.37f);
			DOFShader.SetFloat("u_FocalDepthTemporal", CenterDepthSmooth);
			DOFShader.SetVector2f("u_CameraPlanes", glm::vec2(MainCamera.GetNearPlane(), MainCamera.GetFarPlane()));
			DOFShader.SetBool("u_LargeKernel", LargeKernelDOF);
			DOFShader.SetFloat("u_COCScale", DOFCOCScale);
			DOFShader.SetFloat("u_BlurScale", DOFBlurScale);
			DOFShader.SetFloat("u_CAScale", DOFCAScale);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, DOFFBOInput.GetTexture());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		}

		// Tonemapping ->

		TonemappedFBO.Bind();
		Tonemapper.Use();
		Tonemapper.SetInteger("u_Texture", 0);
		Tonemapper.SetInteger("u_DOFTexture", 1);
		Tonemapper.SetInteger("u_Depth", 2);
		Tonemapper.SetBool("u_FilmGrain", FilmGrainStrength > 0.001f);
		Tonemapper.SetFloat("u_FilmGrainStrength", FilmGrainStrength);
		Tonemapper.SetFloat("u_Time", glfwGetTime());
		Tonemapper.SetFloat("u_Exposure", ComputedExposure* ExposureMultiplier);
		Tonemapper.SetBool("u_HejlBurgess", HejlBurgessTonemap);
		Tonemapper.SetBool("u_DOF", DOF);
		Tonemapper.SetFloat("u_FocalDepthTemporal", CenterDepthSmooth);
		Tonemapper.SetVector2f("u_CameraPlanes", glm::vec2(MainCamera.GetNearPlane(), MainCamera.GetFarPlane()));
		Tonemapper.SetFloat("u_COCScale", DOFCOCScale);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, PostProcessingFBO.GetTexture());
		
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, DOFFBO.GetTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(0));

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		TonemappedFBO.Unbind();

		// ---- CUBE ITEM ---- //


		if (DoSmallItemCube) {

			BlockDatabase::SetTextureArraysFilteringLinear();

			glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);

			TonemappedFBO.Bind();

			CubeItemRenderer.Use();
			CubeItemRenderer.SetFloat("u_Time", glfwGetTime());
			CubeItemRenderer.SetFloat("u_SunVisibility", sun_visibility);
			CubeItemRenderer.SetVector2f("u_Dimensions", glm::vec2(ColoredFBO.GetWidth(), ColoredFBO.GetHeight()));
			CubeItemRenderer.SetVector3f("u_SunDirection", SunDirection);
			CubeItemRenderer.SetVector3f("u_MoonDirection", MoonDirection);
			CubeItemRenderer.SetVector3f("u_StrongerLightDirection", StrongerLightDirection);
			CubeItemRenderer.SetInteger("u_Sky", 0);
			CubeItemRenderer.SetInteger("u_AlbedoTextures", 1);
			CubeItemRenderer.SetInteger("u_NormalTextures", 2);
			CubeItemRenderer.SetInteger("u_PBRTextures", 3);
			CubeItemRenderer.SetInteger("u_EmissiveTextures", 4);
			CubeItemRenderer.SetInteger("u_RandomHDRI", 5);
			CubeItemRenderer.SetInteger("u_RandomHDRIDiffuse", 6);
			CubeItemRenderer.SetInteger("u_BlueNoise", 7);
			CubeItemRenderer.SetInteger("u_HeldBlockID", world->GetCurrentBlock());
			CubeItemRenderer.SetInteger("u_AntialiasLevel", AntialiasItemCubeLevel);
			CubeItemRenderer.SetInteger("u_SpecularSampleBias", ItemCubeSpecularSampleBias);
			CubeItemRenderer.SetBool("u_SimpleLighting", SimpleLightingItemCube);
			CubeItemRenderer.SetInteger("u_GrassBlockProps[0]", VoxelRT::BlockDatabase::GetBlockID("Grass"));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[1]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[2]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[3]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Top));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[4]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[5]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[6]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Front));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[7]", VoxelRT::BlockDatabase::GetBlockTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[8]", VoxelRT::BlockDatabase::GetBlockNormalTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));
			CubeItemRenderer.SetInteger("u_GrassBlockProps[9]", VoxelRT::BlockDatabase::GetBlockPBRTexture("Grass", VoxelRT::BlockDatabase::BlockFaceType::Bottom));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_CUBE_MAP, SkymapMain.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetTextureArray());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetNormalTextureArray());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetPBRTextureArray());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D_ARRAY, VoxelRT::BlockDatabase::GetEmissiveTextureArray());
			glActiveTexture(GL_TEXTURE5);
			glBindTexture(GL_TEXTURE_CUBE_MAP, RandomHDRI.GetID());
			glActiveTexture(GL_TEXTURE6);
			glBindTexture(GL_TEXTURE_CUBE_MAP, RandomHDRIDiffuse.GetID());

			glActiveTexture(GL_TEXTURE7);
			glBindTexture(GL_TEXTURE_2D, BlueNoiseLowResolution.GetTextureID());

			BlockDataStorageBuffer.Bind(0);

			glBindImageTexture(1, TonemappedFBO.GetTexture(), 0, 0, 0, GL_READ_WRITE, GL_RGBA16F);

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			TonemappedFBO.Unbind();

			BlockDatabase::SetTextureArraysFilteringNearest();
		}


		// FXAA 

		bool DoSecondaryFXAA = true && FXAA;

		if (DoSecondaryFXAA) {

			FXAA_Secondary.Use();
			FXAA_FBO.Bind();

			FXAA_Secondary.SetInteger("u_FramebufferTexture", 0);
			FXAA_Secondary.SetInteger("u_PositionTexture", 1);
			FXAA_Secondary.SetInteger("u_NormalTexture", 2);
			FXAA_Secondary.SetInteger("u_BlockIDTex", 3);
			FXAA_Secondary.SetMatrix4("u_InverseView", inv_view);
			FXAA_Secondary.SetMatrix4("u_InverseProjection", inv_projection);
			FXAA_Secondary.SetVector2f("u_Dimensions", glm::vec2(FXAA_FBO.GetWidth(), FXAA_FBO.GetHeight()));
			FXAA_Secondary.SetBool("u_BoostSkyLuminance", StrongerLightDirection == MoonDirection);
			FXAA_Secondary.SetBool("u_ExponentiallyMagnifyColorDifferences", StrongerLightDirection == MoonDirection);
			FXAA_Secondary.SetFloat("u_Time", glfwGetTime());

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, TonemappedFBO.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();

			FXAA_FBO.Unbind();
		}



		// ---- FINAL ----

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		if (ContrastAdaptiveSharpening) {
			FXAA_Final.Bind();
			glViewport(0, 0, app.GetWidth(), app.GetHeight());
		}

		else {
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			glViewport(0, 0, app.GetWidth(), app.GetHeight());
		}

		FinalShader.Use();
		FinalShader.SetInteger("u_FramebufferTexture", 0);
		FinalShader.SetInteger("u_PositionTexture", 1);
		FinalShader.SetInteger("u_NormalTexture", 2);
		FinalShader.SetInteger("u_BlockIDTex", 3);
		FinalShader.SetInteger("u_ColorLUT", 4);
		FinalShader.SetInteger("u_Padding", PIXEL_PADDING);
		FinalShader.SetBool("u_BrutalFXAA", BrutalFXAA);
		FinalShader.SetBool("u_CAS", ContrastAdaptiveSharpening);
		FinalShader.SetMatrix4("u_InverseView", inv_view);
		FinalShader.SetMatrix4("u_InverseProjection", inv_projection);
		FinalShader.SetVector2f("u_Dimensions", glm::vec2(app.GetWidth(), app.GetHeight()));
		FinalShader.SetBool("u_ColorGrading", SelectedColorGradingLUT >= 0);
		FinalShader.SetInteger("u_SelectedLUT", SelectedColorGradingLUT);
		FinalShader.SetBool("u_ColorDither", ColorDither);
		FinalShader.SetBool("FXAA", FXAA);
		FinalShader.SetBool("u_FXAA", FXAA);
		FinalShader.SetBool("u_RenderItemCube", DoSmallItemCube);
		FinalShader.SetBool("u_BoostSkyLuminance", StrongerLightDirection==MoonDirection);
		FinalShader.SetBool("u_ExponentiallyMagnifyColorDifferences", StrongerLightDirection==MoonDirection);
		FinalShader.SetFloat("u_Time", glfwGetTime());
		FinalShader.SetVector4f("u_FocusedOnBlock", HighlightFocusedBlock ? glm::vec4(FocusedOnBlock) : glm::vec4(-1.0f));

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, DoSecondaryFXAA ? FXAA_FBO.GetTexture() : TonemappedFBO.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(3));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(1));

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, InitialTraceFBO->GetTexture(2));
		
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, ColorGradingLUT.GetTextureID());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();

		// CAS

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		if (ContrastAdaptiveSharpening) {
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			glViewport(0, 0, app.GetWidth(), app.GetHeight());

			CAS_Shader.Use();
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			glViewport(0, 0, app.GetWidth(), app.GetHeight());
			CAS_Shader.SetInteger("u_Texture", 0);
			CAS_Shader.SetFloat("u_SharpenAmount", CAS_SharpenAmount);
			CAS_Shader.SetBool("u_ColorGrading", SelectedColorGradingLUT >= 0);
			CAS_Shader.SetBool("u_ColorDither", ColorDither);
			CAS_Shader.SetInteger("u_ColorLUT", 1);
			CAS_Shader.SetInteger("u_SelectedLUT", SelectedColorGradingLUT);
			CAS_Shader.SetBool("FXAA", FXAA);
			CAS_Shader.SetBool("u_FXAA", FXAA);
			CAS_Shader.SetVector4f("u_FocusedOnBlock", HighlightFocusedBlock ? glm::vec4(FocusedOnBlock) : glm::vec4(-1.0f));

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, FXAA_Final.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, ColorGradingLUT.GetTextureID());

			VAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VAO.Unbind();
		} 

		// Particles

		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		if (RenderParticles)
		{
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			world->UpdateParticles(&MainCamera, 
				InitialTraceFBO->GetTexture(0),
				SoftShadows ? (DenoiseSunShadows ? ShadowFiltered.GetTexture() : ShadowTemporalFBO.GetTexture()) : ShadowRawTrace.GetTexture(),
				DiffuseDenoiseFBO.GetTexture(), DiffuseDenoiseFBO.GetTexture(1),
				SunDirection, MainCamera.GetPosition(), 
				glm::vec2(app.GetWidth(), app.GetHeight()), 
				DeltaTime);
			glDisable(GL_BLEND);
		}

		RendererUI.RenderQuad(glm::vec2(floor((float)app.GetWidth() / 2.0f), floor((float)app.GetHeight() / 2.0f)), &Crosshair, &OCamera);

		
		if (MainPlayer.m_Position.y < 2.0f) {
			MainPlayer.m_Position.y = 127.0f;
			MainPlayer.Camera.SetPosition(MainPlayer.m_Position);
		}


		// Prime numbers taken to make sure that there arent any lag spikes.
		// (You want to make sure that too many heavy operations don't take place on the exact same frame) 

		if (app.GetCurrentFrame() % 643 == 0)
		{
			world->GenerateDistanceField();
			world->RebufferLightChunks();
		}

		if (app.GetCurrentFrame() % 60 == 0)
		{
			world->m_ParticleEmitter.CleanUpList();
		}

		if (app.GetCurrentFrame() % 1200 == 0 && PointVolumetricsToggled) {
			std::cout << "\n-- AUTO REUPLOADED VOLUMETRIC VOLUME DATA -- \n";
			Volumetrics::Reupload();
		}


		// Update world
		world->Update(&MainCamera);


		// Finish Frame
		glFinish();
		app.FinishFrame();

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glUseProgram(0);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		std::string title = "Voxel RT | "; title += BlockDatabase::GetBlockName(world->GetCurrentBlock()); title += "     "; 
		title += MainPlayer.Freefly ? "   |  FREEFLY" : "";
		title += MainPlayer.DisableCollisions ? "   |  NO COLLISIONS" : "";
		GLClasses::DisplayFrameRate(app.GetWindow(), title);


		float CurrentTime = glfwGetTime();
		DeltaTime = CurrentTime - Frametime;
		Frametime = glfwGetTime();
		DeltaSum += DeltaTime;
		ModifiedWorld = false;
		AntialiasItemCubeLevel += AntialiasItemCubeLevel % 2;


		// World bounds check
		//glm::vec3 CP = MainCamera.GetPosition(); // Camera position
		//if (MainCamera.GetPosition().y <= 2.0f || CP.x > WORLD_SIZE_X - 2 || CP.x < 2 || CP.z > WORLD_SIZE_Z - 2 || CP.z < 2)
		//{
		//	MainCamera.SetPosition(glm::vec3(WORLD_SIZE_X/2, 75, WORLD_SIZE_Z/2));
		//	MainPlayer.m_Position = MainCamera.GetPosition();
		//	MainPlayer.m_Velocity = glm::vec3(0.0f);
		//	MainPlayer.m_Acceleration = glm::vec3(0.0f);
		//	MainPlayer.m_isOnGround = false;
		//	MainPlayer.InitialCollisionDone = false;
		//	MainPlayer.InitialCollisionDone2 = false;
		//}

		// velocity clamp
		MainPlayer.ClampVelocity();
		
		if (app.GetCurrentFrame() % 10 == 0) {
			MainPlayer.Camera.SetAspect((float)((float)app.GetWidth()) / (float)(app.GetHeight()));
		}
		
		RenderDistance = glm::clamp(RenderDistance, 0, 350);

		// make sure padding is divisible by 2
		if (PIXEL_PADDING % 2 != 0) { PIXEL_PADDING += 1; }

		////std::cout << MainPlayer.InitialCollisionDone;
	}

	SaveWorld(world, world_name);
	SoundManager::Destroy();
	delete world;
	return;
	exit(0);
}

// Forward declaration.
namespace VoxelRT {
	namespace Scope {
		ParticleSystem::ParticleEmitter* GetWorldParticleEmitter() {
			return world?&world->m_ParticleEmitter:nullptr;
		}
	}

}




//////////////////////////
// End of pipeline code.//
//////////////////////////