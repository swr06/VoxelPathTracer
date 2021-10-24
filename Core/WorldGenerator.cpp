#include "WorldGenerator.h"

#include "BlockDatabase.h"



static uint8_t GRASS_ID = 16;
static uint8_t STONE_ID = 32;
static uint8_t DIRT_ID = 48;
static uint8_t SAND_ID = 64;
typedef int Biome;

Biome GetBiome(float chunk_noise)
{
	// Quantize the noise into various levels and frequency

	if (chunk_noise < 130)
	{
		return 0;
	}

	else
	{
		return 1;
	}
}

void SetVerticalBlocks(VoxelRT::World* world, int x, int z, int y_level, int biome)
{
	for (int y = 0; y < y_level; y++)
	{
		if (biome == 1) {
			if (y >= y_level - 1)
			{
				world->SetBlock(x, y, z, { GRASS_ID });
			}

			else if (y >= y_level - 5)
			{
				world->SetBlock(x, y, z, { DIRT_ID });
			}

			else
			{
				world->SetBlock(x, y, z, { STONE_ID });
			}
		}

		else {
			if (y >= y_level - 1)
			{
				world->SetBlock(x, y, z, { SAND_ID });
			}

			else if (y >= y_level - 8)
			{
				world->SetBlock(x, y, z, { SAND_ID });
			}

			else
			{
				world->SetBlock(x, y, z, { STONE_ID });
			}

		}
	}
}



void VoxelRT::GenerateWorld(World* world, bool gen_type)
{
	static FastNoise BiomeGenerator(rand()%10000); // Biome generator (quintized simplex noise)
	static FastNoise NoiseGenerator(rand()%6942); // Simplex fractal

	BiomeGenerator.SetNoiseType(FastNoise::Simplex);

	GRASS_ID = VoxelRT::BlockDatabase::GetBlockID("Grass");
	DIRT_ID = VoxelRT::BlockDatabase::GetBlockID("Dirt");
	STONE_ID = VoxelRT::BlockDatabase::GetBlockID("Stone");
	SAND_ID = VoxelRT::BlockDatabase::GetBlockID("Sand");

	if (gen_type)
	{
		NoiseGenerator.SetNoiseType(FastNoise::SimplexFractal);
		NoiseGenerator.SetFrequency(0.0035);
		NoiseGenerator.SetFractalOctaves(5);

		for (int x = 0; x < WORLD_SIZE_X; x++)
		{
			for (int z = 0; z < WORLD_SIZE_Z; z++)
			{
				float real_x = x;
				float real_z = z;
				float height;
				float h;

				h = (NoiseGenerator.GetNoise(real_x, real_z));
				height = ((h + 1.0f) / 2.0f) * 24.0f;

				float column_noise = BiomeGenerator.GetNoise(real_x/2.0f, real_z/2.0f);
				column_noise = ((column_noise + 1.0f) / 2) * 240;

				SetVerticalBlocks(world, x, z, height+40, GetBiome(column_noise));
			}
		}
	}

	else
	{
		for (int x = 0; x < WORLD_SIZE_X; x++)
		{
			for (int z = 0; z < WORLD_SIZE_Z; z++)
			{
				float real_x = x;
				float real_z = z;

				SetVerticalBlocks(world, x, z, 50,1);
			}
		}
	}
}