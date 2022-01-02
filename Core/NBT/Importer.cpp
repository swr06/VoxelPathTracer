#include "Importer.h"

namespace VoxelRT {
	namespace MCWorldImporter {
		const int HALF_WORLD_X = 384 / 2;
		const int HALF_WORLD_Z = 384 / 2;
		static glm::vec3 ImportOrigin;

		static std::array<uint8_t, 384 * 128 * 384> ReadWorldData;

		bool IsInBounds(const glm::ivec3& Position) {
			if (Position.x <= -HALF_WORLD_X || Position.x >= HALF_WORLD_X || Position.z <= -HALF_WORLD_X || Position.z >= HALF_WORLD_X ||
				Position.y <= 1 || Position.y >= 128) {
				return true;
			}

			return false;
		}

		bool IsInBounds(int x, int y, int z) {
			const glm::ivec3 Position = glm::ivec3(x, y, z);
			if (Position.x <= -HALF_WORLD_X || Position.x >= HALF_WORLD_X || Position.z <= -HALF_WORLD_X || Position.z >= HALF_WORLD_X ||
				Position.y <= 1 || Position.y >= 128) {
				return true;
			}

			return false;
		}

		int ConvertTo1DIDX(int x, int y, int z, const glm::ivec3& Dimensions) {
			return x + y * Dimensions.x + z * Dimensions.x * Dimensions.y;
		}

		int ConvertTo1DIDXWorld(int x, int y, int z) {
			const glm::ivec3 Dimensions = glm::ivec3(WORLD_SIZE_X, WORLD_SIZE_Y, WORLD_SIZE_Z);
			return x + y * Dimensions.x + z * Dimensions.x * Dimensions.y;
		}

		int ConvertTo1DIDX(int x, int y, int z, int SideLength) {
			return x + y * SideLength + z * SideLength * SideLength;
		}



		glm::ivec3 ConvertTo3DIDX(int index, const glm::ivec3& Dimensions) {
			int z = index / (Dimensions.x * Dimensions.y);
			index -= (z * Dimensions.x * Dimensions.y);
			int y = index / Dimensions.x;
			int x = index % Dimensions.x;
			return glm::ivec3(x, y, z);
		}

		glm::ivec3 ConvertTo3DIDX(int index, const int SideLength) {
			int z = index / (SideLength * SideLength);
			index -= (z * SideLength * SideLength);
			int y = index / SideLength;
			int x = index % SideLength;
			return glm::ivec3(x, y, z);
		}

		bool ChunkInBounds(int cx, int cz) {
			bool Valid = cx < -(384 / 32) || cx >(384 / 32) || cz < -(384 / 32) || cz > (384 / 32);
			return !Valid;
		}

		void WriteVoxel(uint8_t voxel, glm::ivec3 Position)
		{
			Position -= glm::ivec3(ImportOrigin);
			Position.x += HALF_WORLD_X;
			Position.z += HALF_WORLD_Z;

			if (Position.y >= 128 || Position.x >= 384 || Position.z >= 384 || Position.y < 0 || Position.x < 0 || Position.z < 0 || voxel == 0) {
				return;
			}

			int index = ConvertTo1DIDXWorld(Position.x, Position.y, Position.z);

			if (index < 0 || index > 384 * 128 * 384) {
				return;
			}

			ReadWorldData[index] = voxel;
		} // 4524 10 937

		void ImportRegionFile(const std::string& Path) {
			FILE* RegionFilePointer = fopen(Path.c_str(), "rb");

			if (!RegionFilePointer)
			{
				throw "Region file not found!";
				return;
			}

			enkiRegionFile RegionFile = enkiRegionFileLoad(RegionFilePointer);

			for (int i = 0; i < ENKI_MI_REGION_CHUNKS_NUMBER; i++)
			{
				enkiNBTDataStream stream;
				enkiInitNBTDataStreamForChunk(RegionFile, i, &stream);

				if (stream.dataLength)
				{
					enkiChunkBlockData aChunk = enkiNBTReadChunk(&stream);
					enkiMICoordinate chunkOriginPos = enkiGetChunkOrigin(&aChunk); // y always 0

					for (int section = 0; section < ENKI_MI_NUM_SECTIONS_PER_CHUNK; ++section)
					{
						if (aChunk.sections[section])
						{
							enkiMICoordinate sectionOrigin = enkiGetChunkSectionOrigin(&aChunk, section);
							glm::ivec3 storeOrigin = glm::ivec3(sectionOrigin.x, sectionOrigin.y, sectionOrigin.z);
							enkiMICoordinate sPos;

							for (sPos.y = 0; sPos.y < ENKI_MI_SIZE_SECTIONS; ++sPos.y)
							{
								for (sPos.z = 0; sPos.z < ENKI_MI_SIZE_SECTIONS; ++sPos.z)
								{
									for (sPos.x = 0; sPos.x < ENKI_MI_SIZE_SECTIONS; ++sPos.x)
									{
										//uint8_t voxel = enkiGetChunkSectionVoxel(&aChunk, section, sPos);
										enkiMIVoxelData ReadVoxel = enkiGetChunkSectionVoxelData(&aChunk, section, sPos);
										glm::ivec3 StoreLoc = storeOrigin + glm::ivec3(sPos.x, sPos.y, sPos.z);
										uint8_t voxel = ReadVoxel.blockID;
										uint8_t dataval = ReadVoxel.dataValue;

										if (dataval == 0) {
											voxel = BlockDatabase::GetIDFromMCID(voxel);
											WriteVoxel(voxel, StoreLoc);
										}
										
									}
								}
							}
						}
					}

					enkiNBTRewind(&stream);

				}

				enkiNBTFreeAllocations(&stream);
			}

			enkiRegionFileFreeAllocations(&RegionFile);

			fclose(RegionFilePointer);

		}

		void ImportWorld(const std::string& Filepath, void* Output, const glm::vec3& origin)
		{
			ImportOrigin = origin;
			memset(&ReadWorldData, 0, 384 * 128 * 384);
			memset(Output, 0, 384 * 128 * 384);

			const std::filesystem::path Directory{ Filepath.c_str() };

			for (auto const& dir_entry : std::filesystem::directory_iterator{ Directory }) {
				std::string File = dir_entry.path().u8string();
				if (File.substr(File.find_last_of(".") + 1) == "mca") {
					ImportRegionFile(File);
				}
			}

			memcpy(Output, &ReadWorldData, 384 * 128 * 384);
		}

	}
}


