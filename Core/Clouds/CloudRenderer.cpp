#include "CloudRenderer.h"

static float CloudResolution = 1.0f;
static bool Checkerboard = true;
static float Coverage = 0.128f;
static float BoxSize = 220.0f;
static float DetailStrength = 0.0f;
static bool Bayer = true;
static bool HighQualityClouds = false;

static std::unique_ptr<GLClasses::Shader> CloudShaderPtr;
static std::unique_ptr<GLClasses::Shader> CloudTemporalFilterPtr;
static std::unique_ptr<GLClasses::Shader> CloudCheckerUpscalerPtr;
static std::unique_ptr<Clouds::NoiseTexture3D> CloudNoise;
static std::unique_ptr<Clouds::NoiseTexture3D> CloudNoiseDetail;

namespace Clouds {
	GLClasses::Framebuffer CurlCloudNoise(128, 128, { { GL_RGB16F, GL_RGB, GL_FLOAT } }, false, false);
}

void Clouds::CloudRenderer::Initialize()
{
	CloudShaderPtr = std::unique_ptr<GLClasses::Shader>(new GLClasses::Shader);
	CloudTemporalFilterPtr = std::unique_ptr<GLClasses::Shader>(new GLClasses::Shader);
	CloudCheckerUpscalerPtr = std::unique_ptr<GLClasses::Shader>(new GLClasses::Shader);
	CloudNoise = std::unique_ptr<Clouds::NoiseTexture3D>(new Clouds::NoiseTexture3D);
	CloudNoiseDetail = std::unique_ptr<Clouds::NoiseTexture3D>(new Clouds::NoiseTexture3D);

	CloudShaderPtr->CreateShaderProgramFromFile("Core/Shaders/Clouds/CloudVert.glsl", "Core/Shaders/Clouds/CloudFrag.glsl");
	CloudShaderPtr->CompileShaders();
	CloudTemporalFilterPtr->CreateShaderProgramFromFile("Core/Shaders/Clouds/FBOVert.glsl", "Core/Shaders/Clouds/TemporalFilter.glsl");
	CloudTemporalFilterPtr->CompileShaders();
	CloudCheckerUpscalerPtr->CreateShaderProgramFromFile("Core/Shaders/Clouds/FBOVert.glsl", "Core/Shaders/Clouds/CheckerUpscaler.glsl");
	CloudCheckerUpscalerPtr->CompileShaders();

	CloudNoise->CreateTexture(160, 160, 160, nullptr); // 32*5
	CloudNoiseDetail->CreateTexture(64, 64, 64, nullptr);

	std::cout << "\nRendering noise textures!\n";
	CurlCloudNoise.CreateFramebuffer();
	CurlCloudNoise.SetSize(128, 128);
	Clouds::RenderNoise(*CloudNoise, 160, false);
	Clouds::RenderNoise(*CloudNoiseDetail, 64, true);
	Clouds::RenderCurlNoise(CurlCloudNoise);

	glBindTexture(GL_TEXTURE_3D, CloudNoise->GetTextureID());
	glGenerateMipmap(GL_TEXTURE_3D);
	glBindTexture(GL_TEXTURE_3D, 0);

	std::cout << "\nRendered noise textures!\n";
}

GLuint Clouds::CloudRenderer::Update(VoxelRT::FPSCamera& MainCamera,
	const glm::mat4& PrevProjection,
	const glm::mat4& PrevView,
	const glm::vec3& CurrentPosition,
	const glm::vec3& PrevPosition,
	GLClasses::VertexArray& VAO,
	const glm::vec3& SunDirection,
	GLuint BlueNoise,
	int AppWidth, int AppHeight, int CurrentFrame, GLuint atmosphere, GLuint pos_tex, glm::vec3 PreviousPosition, GLuint pos_tex_prev, glm::vec2 modifiers, bool Clamp, glm::vec3 DetailParams, float TimeScale)
{
	static CloudFBO CloudFBO_1;
	static CloudFBO CloudFBO_2;
	static GLClasses::Framebuffer CheckerUpscaled(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, false);
	static GLClasses::Framebuffer CloudTemporalFBO1(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, false);
	static GLClasses::Framebuffer CloudTemporalFBO2(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, false);

	glm::mat4 CurrentProjection = MainCamera.GetProjectionMatrix();
	glm::mat4 CurrentView = MainCamera.GetViewMatrix();

	auto& CloudTemporalFBO = (CurrentFrame % 2 == 0) ? CloudTemporalFBO1 : CloudTemporalFBO2;
	auto& PrevCloudTemporalFBO = (CurrentFrame % 2 == 0) ? CloudTemporalFBO2 : CloudTemporalFBO1;
	Clouds::CloudFBO& CloudFBO = (CurrentFrame % 2 == 0) ? CloudFBO_1 : CloudFBO_2;
	Clouds::CloudFBO& PrevCloudFBO = (CurrentFrame % 2 == 0) ? CloudFBO_2 : CloudFBO_1;

	CloudTemporalFBO1.SetSize(AppWidth * 0.75f, AppHeight * 0.75f);
	CloudTemporalFBO2.SetSize(AppWidth * 0.75f, AppHeight * 0.75f);
	float HalfCloudRes = Checkerboard ? (CloudResolution * 0.5f) : CloudResolution;
	CloudFBO.SetDimensions(AppWidth * HalfCloudRes, AppHeight * HalfCloudRes);
	PrevCloudFBO.SetDimensions(AppWidth * HalfCloudRes, AppHeight * HalfCloudRes);
	CheckerUpscaled.SetSize(AppWidth * CloudResolution, AppHeight * CloudResolution);

	GLClasses::Shader& CloudShader = *CloudShaderPtr;
	GLClasses::Shader& TemporalFilter = *CloudTemporalFilterPtr;
	GLClasses::Shader& CheckerUpscaler = *CloudCheckerUpscalerPtr;

	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);

	{
		glm::mat4 inv_view = glm::inverse(MainCamera.GetViewMatrix());
		glm::mat4 inv_projection = glm::inverse(MainCamera.GetProjectionMatrix());

		CloudFBO.Bind();
		CloudShader.Use();

		CloudShader.SetMatrix4("u_InverseView", inv_view);
		CloudShader.SetMatrix4("u_InverseProjection", inv_projection);
		CloudShader.SetInteger("u_WorleyNoise", 0);
		CloudShader.SetInteger("u_CloudNoise", 1);
		CloudShader.SetInteger("u_BlueNoise", 2);
		CloudShader.SetInteger("u_CloudDetailedNoise", 3);
		CloudShader.SetInteger("u_CurlNoise", 6);
		CloudShader.SetFloat("u_Time", glfwGetTime());
		CloudShader.SetFloat("u_Coverage", Coverage);
		CloudShader.SetFloat("BoxSize", BoxSize);
		CloudShader.SetFloat("u_DetailIntensity", DetailStrength);
		CloudShader.SetFloat("u_TimeScale", TimeScale);
		CloudShader.SetInteger("u_CurrentFrame", CurrentFrame);
		CloudShader.SetInteger("u_SliceCount", 256);
		CloudShader.SetVector2f("u_Dimensions", glm::vec2(AppWidth, AppHeight));
		CloudShader.SetVector2f("u_VertDimensions", glm::vec2(AppWidth, AppHeight));
		CloudShader.SetVector2f("u_Modifiers", glm::vec2(modifiers));
		CloudShader.SetVector3f("u_DetailParams", glm::vec3(DetailParams));
		CloudShader.SetVector3f("u_SunDirection", SunDirection);
		CloudShader.SetInteger("u_VertCurrentFrame", CurrentFrame);
		CloudShader.SetInteger("u_Atmosphere", 4);
		CloudShader.SetInteger("u_PositionTex", 5);
		CloudShader.SetBool("u_Checker", Checkerboard);
		CloudShader.SetBool("u_UseBayer", Bayer);
		CloudShader.SetBool("u_HighQualityClouds", HighQualityClouds);
		CloudShader.SetVector2f("u_WindowDimensions", glm::vec2(AppWidth, AppHeight));

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_3D, CloudNoise->GetTextureID());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, BlueNoise);

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_3D, CloudNoiseDetail->GetTextureID());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_CUBE_MAP, atmosphere);

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, pos_tex);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, CurlCloudNoise.GetTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();
	}

	if (Checkerboard)
	{
		CheckerUpscaler.Use();
		CheckerUpscaled.Bind();

		CheckerUpscaler.SetInteger("u_CurrentFrame", CurrentFrame);
		CheckerUpscaler.SetInteger("u_ColorTexture", 0);
		CheckerUpscaler.SetInteger("u_PreviousColorTexture", 2);
		CheckerUpscaler.SetMatrix4("u_PreviousProjection", PrevProjection);
		CheckerUpscaler.SetMatrix4("u_PreviousView", PrevView);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, CloudFBO.GetCloudTexture());
		
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, PrevCloudFBO.GetCloudTexture());

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();
	}

	//// Temporally filter the clouds
	{
		TemporalFilter.Use();
		CloudTemporalFBO.Bind();

		TemporalFilter.SetInteger("u_CurrentColorTexture", 0);
		TemporalFilter.SetInteger("u_PreviousColorTexture", 1);
		TemporalFilter.SetInteger("u_CurrentPositionData", 2);
		TemporalFilter.SetInteger("u_PrevPositionData", 3);
		TemporalFilter.SetMatrix4("u_PrevProjection", PrevProjection);
		TemporalFilter.SetMatrix4("u_PrevView", PrevView);
		TemporalFilter.SetMatrix4("u_Projection", CurrentProjection);
		TemporalFilter.SetMatrix4("u_View", CurrentView);
		TemporalFilter.SetMatrix4("u_InverseProjection", glm::inverse(CurrentProjection));
		TemporalFilter.SetMatrix4("u_InverseView", glm::inverse(CurrentView));
		TemporalFilter.SetFloat("u_Time", glfwGetTime());
		TemporalFilter.SetVector3f("u_CurrentPosition", MainCamera.GetPosition());
		TemporalFilter.SetVector3f("u_PreviousPosition", PreviousPosition);
		TemporalFilter.SetBool("u_Clamp", Clamp);

		float mix_factor = (CurrentPosition != PrevPosition) ? 0.25f : 0.75f;
		TemporalFilter.SetFloat("u_MixModifier", mix_factor);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, Checkerboard ? CheckerUpscaled.GetTexture() : CloudFBO.GetCloudTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, PrevCloudTemporalFBO.GetTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, pos_tex);

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, pos_tex_prev);

		VAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		VAO.Unbind();
	}

	return CloudTemporalFBO.GetTexture();
}

void Clouds::CloudRenderer::RecompileShaders()
{
	CloudShaderPtr->Recompile();
	CloudTemporalFilterPtr->Recompile();
	CloudCheckerUpscalerPtr->Recompile();
}

void Clouds::CloudRenderer::SetChecker(bool v)
{
	Checkerboard = v;
}

void Clouds::CloudRenderer::SetCoverage(float v)
{
	Coverage = v;
}

void Clouds::CloudRenderer::SetBayer(bool v)
{
	Bayer = v;
}

void Clouds::CloudRenderer::SetDetailContribution(float v)
{
	DetailStrength = v;
}

float Clouds::CloudRenderer::GetBoxSize()
{
	return BoxSize;
}

void Clouds::CloudRenderer::SetQuality(bool v)
{
	HighQualityClouds = v;
}

void Clouds::CloudRenderer::SetResolution(float v)
{
	CloudResolution = v;
}
