#include "ShaderManager.h"
#include <sstream>

static std::unordered_map<std::string, GLClasses::Shader> ShaderManager_ShaderMap;

void VoxelRT::ShaderManager::CreateShaders()
{
	AddShader("INITIAL_TRACE", "Core/Shaders/InitialRayTraceVert.glsl", "Core/Shaders/InitialRayTraceFrag.glsl");
	AddShader("FINAL_SHADER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/FBOFrag.glsl");
	AddShader("DIFFUSE_TRACE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/DiffuseRayTraceFrag.glsl");
	AddShader("MAIN_TEMPORAL_FILER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalFilter.glsl");
	AddShader("SHADOW_TRACE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ShadowRayTraceFrag.glsl");
	AddShader("SHADOW_FILTER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ShadowFilter.glsl");
	AddShader("GAUSSIAN_SPATIAL_FILTER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SpatialFilter.glsl");
	AddShader("ATROUS_SPATIAL_FILTER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/AtrousSpatialFilter.glsl");
	
	AddShader("REFLECTION_TRACE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionTraceFrag.glsl");
	AddShader("REFLECTION_DENOISER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionDenoiser.glsl");
	AddShader("REFLECTION_DENOISER_NEW", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionDenoiserNew.glsl");
	AddShader("SPECULAR_TEMPORAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SpecularTemporalFilter.glsl");

	AddShader("RTAO", "Core/Shaders/FBOVert.glsl", "Core/Shaders/RaytracedAO.glsl");
	AddShader("POST_PROCESS", "Core/Shaders/PostProcessingVert.glsl", "Core/Shaders/PostProcessingFrag.glsl");
	AddShader("COLOR_SHADER", "Core/Shaders/ColorPassVert.glsl", "Core/Shaders/ColorPassFrag.glsl");
	AddShader("TEMPORAL_AA", "Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalAA.glsl");

	AddShader("SSAO", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAO.glsl");
	AddShader("BILATERAL_BLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BilateralBlur.glsl");

	AddShader("VARIANCE_ESTIMATOR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SampleVarianceCompute.glsl");
	
	
	AddShader("SIMPLE_DOWNSAMPLE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SimpleDownsampleFrag.glsl");
	AddShader("LUMA_AVERAGER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/CalculateAverageLuminance.glsl");
	AddShader("VOLUMETRIC_SCATTERING", "Core/Shaders/FBOVert.glsl", "Core/Shaders/VolumetricLighting.glsl");
	AddShader("SSAO_BLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAOBlur.glsl");
	AddShader("SPATIAL_INITIAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Spatial3x3Initial.glsl");
	AddShader("CHECKER_RECONSTRUCT", "Core/Shaders/FBOVert.glsl", "Core/Shaders/CheckerboardReconstruct.glsl");
	AddShader("SPECULAR_CHECKER_RECONSTRUCT", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SpecularCheckerReconstruct.glsl");

	// svgf
	AddShader("SVGF_TEMPORAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SVGF/TemporalFilter.glsl");
	AddShader("SVGF_SPATIAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SVGF/SpatialFilter.glsl");
	AddShader("SVGF_VARIANCE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SVGF/VarianceEstimate.glsl");

	// Volumetrics
	AddShader("VOLUMETRICS_COMPUTE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Volumetrics/ComputeVolumetrics.glsl");
	AddShader("VOLUMETRICS_DENOISER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Volumetrics/Denoiser.glsl");
	AddShader("GAUSSIAN_9TAP_OPTIMIZED", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Gaussian9TapSinglePass.glsl");
	AddShader("GAUSSIAN_5TAP_OPTIMIZED", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Gaussian5TapSinglePass.glsl");
	
	// CAS
	AddShader("CONTRAST_ADAPTIVE_SHARPENING", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ContrastAdaptiveSharpening.glsl");

	// FXAA
	AddShader("FXAA_SECONDARY", "Core/Shaders/FBOVert.glsl", "Core/Shaders/FXAA311.glsl");

	AddShader("DOF", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BokehDOF.glsl");
	AddShader("BICUBIC_DOWNSAMPLE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BicubicDownsample.glsl");

	// anti flicker
	AddShader("ANTI_FLICKER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/HandleFlicker.glsl");
	AddShader("BILATERAL_HITDIST1", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BlurHitDistance_Bilateral.glsl");
	AddShader("BILATERAL_HITDIST2", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BlurHitDistancePass2_Bilateral.glsl");
	AddShader("SSSBLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/FakeSSSSS.glsl");
	AddShader("SSSBLUR_POISSON", "Core/Shaders/FBOVert.glsl", "Core/Shaders/FakeSSSSS_Poisson.glsl");
	AddShader("CUBE_ITEM_RENDERER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/CubeItemShader.glsl");
	AddShader("DEPTH_DOWNSAMPLE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/DownsampleDepth4x.glsl");
	AddShader("NORMAL_DOWNSAMPLE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/DownsampleNormals4x.glsl");
	AddShader("GBUFFER_DOWNSAMPLE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/GBufferDownsampler.glsl");
	AddShader("GBUFFER_GENERATE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/GenerateGBuffer.glsl");
	AddShader("TONEMAPPER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Tonemapper.glsl");
	AddShader("BLOOM_COMBINE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/CombineBloom.glsl");
	AddShader("DIFFRACTION_SPIKES", "Core/Shaders/FBOVert.glsl", "Core/Shaders/DiffractionSpikes.glsl");


}

void VoxelRT::ShaderManager::AddShader(const std::string& name, const std::string& vert, const std::string& frag, const std::string& geo)
{
	auto exists = ShaderManager_ShaderMap.find(name);

	if (exists == ShaderManager_ShaderMap.end())
	{
		ShaderManager_ShaderMap.emplace(name, GLClasses::Shader());
		ShaderManager_ShaderMap.at(name).CreateShaderProgramFromFile(vert, frag);
		ShaderManager_ShaderMap.at(name).CompileShaders();
	}

	else
	{
		std::string err = "A shader with the name : " + name + "  already exists!";
		throw err;
	}
}

GLClasses::Shader& VoxelRT::ShaderManager::GetShader(const std::string& name)
{
	auto exists = ShaderManager_ShaderMap.find(name);

	if (exists != ShaderManager_ShaderMap.end())
	{
		return ShaderManager_ShaderMap.at(name);
	}

	else
	{
		throw "Shader that doesn't exist trying to be accessed!";
	}
}

GLuint VoxelRT::ShaderManager::GetShaderID(const std::string& name)
{
	auto exists = ShaderManager_ShaderMap.find(name);

	if (exists != ShaderManager_ShaderMap.end())
	{
		return ShaderManager_ShaderMap.at(name).GetProgramID();
	}

	else
	{
		throw "Shader that doesn't exist trying to be accessed!";
	}
}

void VoxelRT::ShaderManager::RecompileShaders()
{
	for (auto& e : ShaderManager_ShaderMap)
	{
		e.second.Recompile();
	}
}
