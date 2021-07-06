#include "BlockDataSSBO.h"

namespace VoxelRT
{
	void BlockDataSSBO::CreateBuffers()
	{
		m_SSBO = 0;

		std::array<int, 128> AlbedoData;
		std::array<int, 128> NormalData;
		std::array<int, 128> PBRData;
		std::array<int, 128> EmissiveData;
		std::array<int, 128> Transparent;
		std::vector<int> TotalData;

		for (int b = 0; b < 128; b++)
		{
			uint8_t i = static_cast<uint8_t>(b);
			AlbedoData[i] = BlockDatabase::GetBlockTexture({ i }, BlockDatabase::BlockFaceType::Front);
			NormalData[i] = BlockDatabase::GetBlockNormalTexture({ i }, BlockDatabase::BlockFaceType::Front);
			PBRData[i] = BlockDatabase::GetBlockPBRTexture({ i }, BlockDatabase::BlockFaceType::Front);
			EmissiveData[i] = BlockDatabase::GetBlockEmissiveTexture({ i });
			Transparent[i] = BlockDatabase::IsBlockTransparent({ i }) ? 1 : 0;
		}
		
		int TotalSize = (128 * sizeof(int)) * 5;

		TotalData.insert(TotalData.end(), std::begin(AlbedoData), std::end(AlbedoData));
		TotalData.insert(TotalData.end(), std::begin(NormalData), std::end(NormalData));
		TotalData.insert(TotalData.end(), std::begin(PBRData), std::end(PBRData));
		TotalData.insert(TotalData.end(), std::begin(EmissiveData), std::end(EmissiveData));
		TotalData.insert(TotalData.end(), std::begin(Transparent), std::end(Transparent));

		glGenBuffers(1, &m_SSBO);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_SSBO);
		glBufferData(GL_SHADER_STORAGE_BUFFER, TotalSize, TotalData.data(), GL_STATIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}

	void BlockDataSSBO::Bind(int idx)
	{
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, idx, m_SSBO);
	}
}