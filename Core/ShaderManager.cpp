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
	AddShader("SPATIAL_FILTER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SpatialFilter.glsl");
	
	AddShader("REFLECTION_TRACE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionTraceFrag.glsl");
	AddShader("REFLECTION_DENOISER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/ReflectionDenoiser.glsl");
	AddShader("SPECULAR_TEMPORAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SpecularTemporalFilter.glsl");

	AddShader("RTAO", "Core/Shaders/FBOVert.glsl", "Core/Shaders/RaytracedAO.glsl");
	AddShader("POST_PROCESS", "Core/Shaders/PostProcessingVert.glsl", "Core/Shaders/PostProcessingFrag.glsl");
	AddShader("COLOR_SHADER", "Core/Shaders/ColorPassVert.glsl", "Core/Shaders/ColorPassFrag.glsl");
	AddShader("TEMPORAL_AA", "Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalAA.glsl");

	AddShader("SSAO", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAO.glsl");
	AddShader("BILATERAL_BLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BilateralBlur.glsl");


	AddShader("SMART_DENOISER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SmartDenoise.glsl");
	AddShader("SIMPLE_DOWNSAMPLE", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SimpleDownsampleFrag.glsl");
	AddShader("LUMA_AVERAGER", "Core/Shaders/FBOVert.glsl", "Core/Shaders/CalculateAverageLuminance.glsl");
	AddShader("VOLUMETRIC_SCATTERING", "Core/Shaders/FBOVert.glsl", "Core/Shaders/VolumetricLighting.glsl");
	AddShader("SSAO_BLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAOBlur.glsl");
	AddShader("SPATIAL_INITIAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Spatial2x2Initial.glsl");
	AddShader("CHECKER_RECONSTRUCT", "Core/Shaders/FBOVert.glsl", "Core/Shaders/CheckerboardReconstruct.glsl");
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
