#include "ShaderManager.h"
#include <sstream>

static std::unordered_map<std::string, GLClasses::Shader> ShaderManager_ShaderMap;

void Blocks::ShaderManager::CreateShaders()
{
	AddShader("RENDER_SHADER", "Core/Shaders/BlockVert.glsl", "Core/Shaders/BlockFrag.glsl");
	AddShader("POST_PROCESSING", "Core/Shaders/Tonemapping/TonemapVert.glsl", "Core/Shaders/Tonemapping/ACES.glsl");
	AddShader("VOLUMETRICS", "Core/Shaders/VolumetricLightingVert.glsl", "Core/Shaders/VolumetricLightingFrag.glsl");
	AddShader("SSR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSRFrag.glsl");
	AddShader("SS_REFRACTIONS", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSRefractionsFrag.glsl");
	AddShader("WATER", "Core/Shaders/WaterVert.glsl", "Core/Shaders/WaterFrag.glsl");
	AddShader("DEPTH_PREPASS", "Core/Shaders/DepthPrepassVert.glsl", "Core/Shaders/DepthPrepassFrag.glsl");
	AddShader("ATMOSPHERE_COMBINE", "Core/Shaders/AtmosphereCombineVert.glsl", "Core/Shaders/AtmosphereCombineFrag.glsl");
	AddShader("ATMOSPHERE", "Core/Shaders/AtmosphereVertex.glsl", "Core/Shaders/AtmosphereFrag.glsl");
	AddShader("BILATERAL_BLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BilateralBlurFrag.glsl");
	AddShader("GAUSSIAN_VERTICAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/GaussianBlurVertical.glsl");
	AddShader("GAUSSIAN_HORIZONTAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/GaussianBlurHorizontal.glsl");
	AddShader("BLOOM_BRIGHT", "Core/Shaders/FBOVert.glsl", "Core/Shaders/BloomMaskFrag.glsl");
	AddShader("CUBEMAP_REFLECTION", "Core/Shaders/CubemapReflectionVert.glsl", "Core/Shaders/CubemapReflectionFrag.glsl");
	AddShader("DEPTH", "Core/Shaders/DepthVert.glsl", "Core/Shaders/DepthFrag.glsl");
	AddShader("SSAO", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAO.glsl");
	AddShader("SSAO_BLUR", "Core/Shaders/FBOVert.glsl", "Core/Shaders/SSAOBlur.glsl");
	AddShader("GAUSSIAN_SINGLEPASS", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Gaussian9TapSinglePass.glsl");
	AddShader("GAUSSIAN_SINGLEPASS_5TAP", "Core/Shaders/FBOVert.glsl", "Core/Shaders/Gaussian5TapSinglePass.glsl");
	AddShader("TAA", "Core/Shaders/FBOVert.glsl", "Core/Shaders/TemporalAA.glsl");
	AddShader("FINAL", "Core/Shaders/FBOVert.glsl", "Core/Shaders/FinalFrag.glsl");
}

void Blocks::ShaderManager::AddShader(const std::string& name, const std::string& vert, const std::string& frag, const std::string& geo)
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

GLClasses::Shader& Blocks::ShaderManager::GetShader(const std::string& name)
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

GLuint Blocks::ShaderManager::GetShaderID(const std::string& name)
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

void Blocks::ShaderManager::RecompileShaders()
{
	for (auto& e : ShaderManager_ShaderMap)
	{
		e.second.Recompile();
	}
}
