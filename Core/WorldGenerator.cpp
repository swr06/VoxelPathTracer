#include "WorldGenerator.h"

#include "BlockDatabase.h"

static FastNoise NoiseGenerator(2384); // World generator

uint8_t GRASS_ID = 16;
uint8_t STONE_ID = 32;
uint8_t DIRT_ID = 48;

void SetVerticalBlocks(VoxelRT::World* world, int x, int z, int y_level)
{
	for (int y = 0; y < y_level; y++)
	{
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
}


void VoxelRT::GenerateWorld(World* world, bool gen_type)
{
	GRASS_ID = VoxelRT::BlockDatabase::GetBlockID("Grass");
	DIRT_ID = VoxelRT::BlockDatabase::GetBlockID("Dirt");
	STONE_ID = VoxelRT::BlockDatabase::GetBlockID("Stone");

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
				height = ((h + 1.0f) / 2.0f) * 100;

				SetVerticalBlocks(world, x, z, height);
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

				SetVerticalBlocks(world, x, z, 40);
			}
		}
	}
}