#include "VolumetricFloodFill.h"

#include <memory>

#define PACK_U16(lsb, msb) ((uint16_t) ( ((uint16_t)(lsb) & 0xFF) | (((uint16_t)(msb) & 0xFF) << 8) ))

// Used to debug flood fill algorithm
//#define VOXEL_RT_VOLUMETRICS_DEBUG

namespace VoxelRT {
	
	static GLuint VolumetricVolume = 0;
	static GLuint AverageColorSSBO = 0;
	static std::queue<LightNode> LightBFS;
	static std::queue<LightRemovalNode> LightRemovalBFS;
	static World* VolumetricWorldPtr = nullptr;
	static std::unique_ptr<std::array<uint16_t, 384 * 128 * 384>> LightWorldData;

	bool InVoxelVolume(const glm::ivec3& x) {

		if (x.x > 0 && x.y > 0 && x.z > 0 && x.x < 384 && x.y < 128 && x.z < 384)
		{
			return true;
		}

		return false;
	}
};

void VoxelRT::Volumetrics::CreateVolume(World* world, GLuint SSBO_Blockdata, GLuint AlbedoArray)
{
	VolumetricVolume = 0;
	AverageColorSSBO = 0;
	VolumetricWorldPtr = world;

	glGenTextures(1, &VolumetricVolume);
	glBindTexture(GL_TEXTURE_3D, VolumetricVolume);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_R16UI, 384, 128, 384, 0, GL_RED_INTEGER, GL_UNSIGNED_SHORT, nullptr);

	// clear (doesnt automatically clear on some gpus)
	GLClasses::ComputeShader ClearShader;
	ClearShader.CreateComputeShader("Core/Shaders/Volumetrics/ClearData.comp");
	ClearShader.Compile();
	ClearShader.Use();
	glBindImageTexture(0, VolumetricVolume, 0, GL_TRUE, 0, GL_READ_WRITE, GL_R16UI);
	glDispatchCompute(384 / 8, 128 / 8, 384 / 8);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	// Create light data array
	LightWorldData = std::unique_ptr<std::array<uint16_t, 384 * 128 * 384>>(new std::array<uint16_t, 384 * 128 * 384>);

	// memset 
	memset(LightWorldData.get()->data(), 0, 384 * 128 * 384 * sizeof(uint16_t));

	// initialize data ssbo 
	glGenBuffers(1, &AverageColorSSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, AverageColorSSBO);
	int TotalBufferSize = (sizeof(GLfloat) * 4) * 128;
	glBufferData(GL_SHADER_STORAGE_BUFFER, TotalBufferSize, nullptr, GL_STATIC_DRAW);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	// Write data
	GLClasses::ComputeShader DataShader;
	DataShader.CreateComputeShader("Core/Shaders/Volumetrics/PrecomputeAverageBlockColor.comp");
	DataShader.Compile();

	// Dispatch!
	DataShader.Use();
	DataShader.SetInteger("u_BlockAlbedo", 4);
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D_ARRAY, AlbedoArray);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, SSBO_Blockdata);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, AverageColorSSBO);
	glDispatchCompute(1, 1, 1);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
	glFinish();
}

uint8_t VoxelRT::Volumetrics::GetLightValue(const glm::ivec3& p)
{
#ifdef VOXEL_RT_VOLUMETRICS_DEBUG
	std::cout << "GetLightValue() Called";
#endif
	if (!VoxelRT::InVoxelVolume(p)) {
		throw "wtf";
	}

	int idx = p.x + p.y * WORLD_SIZE_X + p.z * WORLD_SIZE_X * WORLD_SIZE_Y;
	auto& arr = *LightWorldData;
	return (arr.at(idx)) & 0xFF;
}

uint8_t VoxelRT::Volumetrics::GetBlockTypeLightValue(const glm::ivec3& p)
{
#ifdef VOXEL_RT_VOLUMETRICS_DEBUG
	std::cout << "GetBlockTypeLightValue() Called";
#endif

	if (!VoxelRT::InVoxelVolume(p)) {
		throw "wtf";
	}

	int idx = p.x + p.y * WORLD_SIZE_X + p.z * WORLD_SIZE_X * WORLD_SIZE_Y;
	auto& arr = *LightWorldData;
	return ((arr.at(idx)) >> 8) & 0xFF ;
}

void VoxelRT::Volumetrics::SetLightValue(const glm::ivec3& p, uint8_t v, uint8_t block)
{
#ifdef VOXEL_RT_VOLUMETRICS_DEBUG
	std::cout << "SetLightValue() Called";
#endif

	if (!VoxelRT::InVoxelVolume(p)) {
		throw "wtf";
	}

	auto& arr = *LightWorldData;
	int idx = p.x + p.y * WORLD_SIZE_X + p.z * WORLD_SIZE_X * WORLD_SIZE_Y;
	arr[idx] = PACK_U16(v, block);
}

void VoxelRT::Volumetrics::UploadLight(const glm::ivec3& p, uint8_t v, uint8_t block, bool should_bind)
{
#ifdef VOXEL_RT_VOLUMETRICS_DEBUG
	std::cout << "UploadLight() Called";
#endif

	if (!VoxelRT::InVoxelVolume(p)) {
		throw "wtf";
	}

	uint16_t data = PACK_U16(v, block);

	if (should_bind) {
		glBindTexture(GL_TEXTURE_3D, VolumetricVolume);
	}

	glTexSubImage3D(GL_TEXTURE_3D, 0, (int)p.x, (int)p.y, (int)p.z, 1, 1, 1, GL_RED_INTEGER, GL_UNSIGNED_SHORT, &data);
}

void VoxelRT::Volumetrics::AddLightToVolume(const glm::ivec3& p, uint8_t block)
{
	VoxelRT::Volumetrics::SetLightValue(glm::ivec3(
		floor(p.x),
		floor(p.y),
		floor(p.z)), 7, block); 
	VoxelRT::Volumetrics::UploadLight(glm::ivec3(
		floor(p.x),
		floor(p.y),
		floor(p.z)), 7, block, true);
	LightBFS.push(LightNode(glm::vec3(
		floor(p.x),
		floor(p.y),
		floor(p.z))));
}

void VoxelRT::Volumetrics::Reupload()
{
	glBindTexture(GL_TEXTURE_3D, VolumetricVolume);
	glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, 384, 128, 384, GL_RED_INTEGER, GL_UNSIGNED_SHORT, LightWorldData.get()->data());
	glBindTexture(GL_TEXTURE_3D, 0);
}

GLuint VoxelRT::Volumetrics::GetVolume()
{
	return VolumetricVolume;
}

GLuint VoxelRT::Volumetrics::GetAverageColorSSBO()
{
	return AverageColorSSBO;
}

std::queue<VoxelRT::LightNode>& VoxelRT::Volumetrics::GetLightBFSQueue()
{
	return LightBFS;
}

std::queue<VoxelRT::LightRemovalNode>& VoxelRT::Volumetrics::GetLightRemovalBFSQueue()
{
	return LightRemovalBFS;
}

void VoxelRT::Volumetrics::PropogateVolume()
{
	World* world = VolumetricWorldPtr;

	glBindTexture(GL_TEXTURE_3D, VolumetricVolume);

	while (!LightBFS.empty())
	{
		LightNode node = VoxelRT::LightBFS.front();
		LightBFS.pop();

		glm::ivec3 pos = node.m_Position;
		uint8_t current_light = VoxelRT::Volumetrics::GetLightValue(pos);
		uint8_t current_block_type = VoxelRT::Volumetrics::GetBlockTypeLightValue(pos);
		glm::ivec3 temp_pos = glm::vec3(0.0f);

		if (temp_pos = glm::vec3(pos.x + 1, pos.y, pos.z); world->GetBlock(temp_pos).block == 0 && Volumetrics::GetLightValue(temp_pos) + 2 < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, current_light - 1, current_block_type);
			Volumetrics::UploadLight(temp_pos, current_light - 1, current_block_type, false);
			LightBFS.push(LightNode(temp_pos));
		}

		if (temp_pos = glm::vec3(pos.x - 1, pos.y, pos.z); world->GetBlock(temp_pos).block == 0 && Volumetrics::GetLightValue(temp_pos) + 2 < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, current_light - 1, current_block_type);
			Volumetrics::UploadLight(temp_pos, current_light - 1, current_block_type, false);
			LightBFS.push(LightNode(temp_pos));
		}

		if (temp_pos = glm::vec3(pos.x, pos.y + 1, pos.z); world->GetBlock(temp_pos).block == 0 && Volumetrics::GetLightValue(temp_pos) + 2 < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, current_light - 1, current_block_type);
			Volumetrics::UploadLight(temp_pos, current_light - 1, current_block_type, false);
			LightBFS.push(LightNode(temp_pos));
		}

		if (temp_pos = glm::vec3(pos.x, pos.y - 1, pos.z); world->GetBlock(temp_pos).block == 0 && Volumetrics::GetLightValue(temp_pos) + 2 < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, current_light - 1, current_block_type);
			Volumetrics::UploadLight(temp_pos, current_light - 1, current_block_type, false);
			LightBFS.push(LightNode(temp_pos));
		}

		if (temp_pos = glm::vec3(pos.x, pos.y, pos.z - 1); world->GetBlock(temp_pos).block == 0 && Volumetrics::GetLightValue(temp_pos) + 2 < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, current_light - 1, current_block_type);
			Volumetrics::UploadLight(temp_pos, current_light - 1, current_block_type, false);
			LightBFS.push(LightNode(temp_pos));
		}

		if (temp_pos = glm::vec3(pos.x, pos.y, pos.z + 1); world->GetBlock(temp_pos).block == 0 && Volumetrics::GetLightValue(temp_pos) + 2 < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, current_light - 1, current_block_type);
			Volumetrics::UploadLight(temp_pos, current_light - 1, current_block_type, false);
			LightBFS.push(LightNode(temp_pos));
		}
	}
}

void VoxelRT::Volumetrics::DepropogateVolume()
{
	glBindTexture(GL_TEXTURE_3D, VolumetricVolume);

	while (!LightRemovalBFS.empty())
	{
		LightRemovalNode node = LightRemovalBFS.front();
		LightRemovalBFS.pop();

		glm::ivec3 pos = node.m_Position;
		uint8_t current_light = node.m_LightValue;
		glm::ivec3 temp_pos = glm::vec3(0.0f);
		uint8_t neighbouring_light;

		// x + 1
		temp_pos = glm::vec3(pos.x + 1, pos.y, pos.z);
		neighbouring_light = VoxelRT::Volumetrics::GetLightValue(temp_pos);
		uint8_t CurrentBlock = VoxelRT::Volumetrics::GetBlockTypeLightValue(temp_pos);

		if (neighbouring_light != 0 && neighbouring_light < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, 0, CurrentBlock);
			Volumetrics::UploadLight(temp_pos, 0, CurrentBlock, false);
			LightRemovalBFS.push(LightRemovalNode(temp_pos, neighbouring_light));
		}

		else if (neighbouring_light >= current_light)
		{
			LightBFS.push(LightNode(temp_pos));
		}

		// x - 1

		temp_pos = glm::vec3(pos.x - 1, pos.y, pos.z);
		neighbouring_light = VoxelRT::Volumetrics::GetLightValue(temp_pos);
		CurrentBlock = VoxelRT::Volumetrics::GetBlockTypeLightValue(temp_pos);

		if (neighbouring_light != 0 && neighbouring_light < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, 0, CurrentBlock);
			Volumetrics::UploadLight(temp_pos, 0, CurrentBlock, false);
			LightRemovalBFS.push(LightRemovalNode(temp_pos, neighbouring_light));
		}

		else if (neighbouring_light >= current_light)
		{
			LightBFS.push(LightNode(temp_pos));
		}

		// y + 1

		temp_pos = glm::vec3(pos.x, pos.y + 1, pos.z);
		neighbouring_light = VoxelRT::Volumetrics::GetLightValue(temp_pos);
		CurrentBlock = VoxelRT::Volumetrics::GetBlockTypeLightValue(temp_pos);

		if (neighbouring_light != 0 && neighbouring_light < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, 0, CurrentBlock);
			Volumetrics::UploadLight(temp_pos, 0, CurrentBlock, false);
			LightRemovalBFS.push(LightRemovalNode(temp_pos, neighbouring_light));
		}

		else if (neighbouring_light >= current_light)
		{
			LightBFS.push(LightNode(temp_pos));
		}

		// y - 1

		temp_pos = glm::vec3(pos.x, pos.y - 1, pos.z);
		neighbouring_light = VoxelRT::Volumetrics::GetLightValue(temp_pos);
		CurrentBlock = VoxelRT::Volumetrics::GetBlockTypeLightValue(temp_pos);

		if (neighbouring_light != 0 && neighbouring_light < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, 0, CurrentBlock);
			Volumetrics::UploadLight(temp_pos, 0, CurrentBlock, false);
			LightRemovalBFS.push(LightRemovalNode(temp_pos, neighbouring_light));
		}

		else if (neighbouring_light >= current_light)
		{
			LightBFS.push(LightNode(temp_pos));
		}

		// z + 1

		temp_pos = glm::vec3(pos.x, pos.y, pos.z + 1);
		neighbouring_light = VoxelRT::Volumetrics::GetLightValue(temp_pos);
		CurrentBlock = VoxelRT::Volumetrics::GetBlockTypeLightValue(temp_pos);

		if (neighbouring_light != 0 && neighbouring_light < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, 0, CurrentBlock);
			Volumetrics::UploadLight(temp_pos, 0, CurrentBlock, false);
			LightRemovalBFS.push(LightRemovalNode(temp_pos, neighbouring_light));
		}

		else if (neighbouring_light >= current_light)
		{
			LightBFS.push(LightNode(temp_pos));
		}

		// z - 1

		temp_pos = glm::vec3(pos.x, pos.y, pos.z - 1);
		neighbouring_light = VoxelRT::Volumetrics::GetLightValue(temp_pos);
		CurrentBlock = VoxelRT::Volumetrics::GetBlockTypeLightValue(temp_pos);

		if (neighbouring_light != 0 && neighbouring_light < current_light)
		{
			Volumetrics::SetLightValue(temp_pos, 0, CurrentBlock);
			Volumetrics::UploadLight(temp_pos, 0, CurrentBlock, false);
			LightRemovalBFS.push(LightRemovalNode(temp_pos, neighbouring_light));
		}

		else if (neighbouring_light >= current_light)
		{
			LightBFS.push(LightNode(temp_pos));
		}
	}
}

