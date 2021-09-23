#include "SoundManager.h"

#include "Audio.h"

#include <glfw/glfw3.h>

extern float VoxelRT_VolumeMultiplier;

namespace VoxelRT
{
	static irrklang::ISoundEngine* MainSoundEngine;
	static std::unordered_map<std::string, std::vector<std::string>> SoundPathsStep;
	static std::unordered_map<std::string, std::vector<std::string>> SoundPathsModify;
	static double PreviousStepTime = -1.0f;

	void SoundManager::InitializeSoundManager()
	{
		MainSoundEngine = irrklang::createIrrKlangDevice();

		if (!MainSoundEngine)
		{
			throw "Sound engine couldnt be initialized :/";
		}

		SetPack(true);
	}

	void SoundManager::UpdatePosition(const glm::vec3& Front, const glm::vec3& Position, const glm::vec3& Up)
	{
		const glm::vec3& position = Position;
		const glm::vec3& front = Front;
		const glm::vec3& up = Up;
		irrklang::vec3df _position(position.x, position.y, position.z);
		irrklang::vec3df _vps(0, 0, 0); // not really needed
		irrklang::vec3df _front(front.x, front.y, front.z);
		irrklang::vec3df _up(up.x, up.y, up.z);
		MainSoundEngine->setListenerPosition(_position, _front, _vps, _up);
	}

	void SoundManager::PlaySound(const std::string& s, const glm::vec3& p, float d, float v, bool pause)
	{
		Audio::Audio3D aud(s, p, MainSoundEngine, false);
		aud.p_Sound->setMinDistance(d);
		aud.p_Sound->setVolume(v);
		aud.p_Sound->setIsPaused(pause);
	}

	void SoundManager::PlayBlockSound(uint8_t block, const glm::vec3& p, bool type)
	{
		// Air.
		if (block <= 0)
		{
			return;
		}

		if (PreviousStepTime >= 0.0f)
		{
			float DeltaStep = glfwGetTime() - PreviousStepTime;

			if (DeltaStep < 0.3750f && type) {
				return;
			}
		}

		// 1 : step, 0 : modify
		std::string Type = type ? BlockDatabase::GetStepSound(block) : BlockDatabase::GetModifySound(block);
		auto& list = type ? SoundPathsStep : SoundPathsModify;

		if (list.find(Type) == list.end()) {
			std::cout << "\nCant play sound for block : " << block << "\n" ;
			return;
		}

		int rand_idx = rand() % list[Type].size();
		rand_idx = glm::clamp(rand_idx, 0, 8);
		std::string Path = list[Type][rand_idx];

		float v = type ? 6.0f : 6.0f; v *= VoxelRT_VolumeMultiplier;
		PlaySound(Path, p, 0.8675f, v, false);

		PreviousStepTime = type ? glfwGetTime() : PreviousStepTime;
	}

	void SoundManager::Destroy()
	{
		if (!MainSoundEngine) 
		{ 
			return; 
		}

		MainSoundEngine->drop();
	}

	void SoundManager::SetPack(bool TYPE)
	{
		if (TYPE) {
			// Step :
			std::string BaseStepPath = "Res/Sounds/Default/Step/";
			SoundPathsStep["Cloth"] = { BaseStepPath + "cloth1.ogg", BaseStepPath + "cloth2.ogg", BaseStepPath + "cloth3.ogg", BaseStepPath + "cloth4.ogg" };

			SoundPathsStep["Sand"] = { BaseStepPath + "sand1.ogg", BaseStepPath + "sand2.ogg", BaseStepPath + "sand3.ogg", BaseStepPath + "sand4.ogg" };

			SoundPathsStep["Snow"] = { BaseStepPath + "snow_walk1.ogg", BaseStepPath + "snow_walk2.ogg", BaseStepPath + "snow_walk3.ogg", BaseStepPath + "snow_walk4.ogg",
										BaseStepPath + "snow_walk5.ogg", BaseStepPath + "snow_walk6.ogg", BaseStepPath + "snow_walk7.ogg", BaseStepPath + "snow_walk8.ogg" };

			SoundPathsStep["Grass"] = { BaseStepPath + "grass_walk1.ogg", BaseStepPath + "grass_walk2.ogg", BaseStepPath + "grass_walk3.ogg", BaseStepPath + "grass_walk4.ogg",
										BaseStepPath + "grass_walk5.ogg", BaseStepPath + "grass_walk6.ogg", BaseStepPath + "grass_walk7.ogg", BaseStepPath + "grass_walk8.ogg" };

			SoundPathsStep["Dirt"] = { BaseStepPath + "dirt_walk1.ogg", BaseStepPath + "dirt_walk2.ogg", BaseStepPath + "dirt_walk3.ogg", BaseStepPath + "dirt_walk4.ogg",
										BaseStepPath + "dirt_walk5.ogg", BaseStepPath + "dirt_walk6.ogg", BaseStepPath + "dirt_walk7.ogg", BaseStepPath + "dirt_walk8.ogg" };

			SoundPathsStep["Leaves"] = { BaseStepPath + "leaves_through1.ogg", BaseStepPath + "leaves_through2.ogg", BaseStepPath + "leaves_through3.ogg", BaseStepPath + "leaves_through4.ogg",
										BaseStepPath + "leaves_through5.ogg", BaseStepPath + "leaves_through6.ogg", BaseStepPath + "leaves_through7.ogg" };

			SoundPathsStep["Gravel"] = { BaseStepPath + "gravel_walk1.ogg", BaseStepPath + "gravel_walk2.ogg", BaseStepPath + "gravel_walk3.ogg", BaseStepPath + "gravel_walk4.ogg",
										 BaseStepPath + "gravel_walk5.ogg", BaseStepPath + "gravel_walk6.ogg", BaseStepPath + "gravel_walk7.ogg", BaseStepPath + "gravel_walk8.ogg" };


			SoundPathsStep["Stone"] = { BaseStepPath + "stone_walk1.ogg", BaseStepPath + "stone_walk2.ogg", BaseStepPath + "stone_walk3.ogg", BaseStepPath + "stone_walk4.ogg",
										BaseStepPath + "stone_walk5.ogg", BaseStepPath + "stone_walk6.ogg", BaseStepPath + "stone_walk7.ogg", BaseStepPath + "stone_walk8.ogg" };

			SoundPathsStep["Wood"] = { BaseStepPath + "wood_walk1.ogg", BaseStepPath + "wood_walk2.ogg", BaseStepPath + "wood_walk3.ogg", BaseStepPath + "wood_walk4.ogg",
									   BaseStepPath + "wood_walk5.ogg", BaseStepPath + "wood_walk6.ogg", BaseStepPath + "wood_walk7.ogg", BaseStepPath + "wood_walk8.ogg" };

			SoundPathsStep["Metal"] = { BaseStepPath + "metalbar_walk1.ogg", BaseStepPath + "metalbar_walk2.ogg", BaseStepPath + "metalbar_walk3.ogg", BaseStepPath + "metalbar_walk4.ogg",
										BaseStepPath + "metalbar_walk5.ogg", BaseStepPath + "metalbar_walk6.ogg", BaseStepPath + "metalbar_walk7.ogg", BaseStepPath + "metalbar_walk8.ogg" };

			SoundPathsStep["Marble"] = { BaseStepPath + "marble_walk1.ogg", BaseStepPath + "marble_walk2.ogg", BaseStepPath + "marble_walk3.ogg", BaseStepPath + "marble_walk4.ogg",
										 BaseStepPath + "marble_walk5.ogg", BaseStepPath + "marble_walk6.ogg", BaseStepPath + "marble_walk7.ogg", BaseStepPath + "marble_walk8.ogg" };

			SoundPathsStep["Concrete"] = { BaseStepPath + "concrete_walk1.ogg", BaseStepPath + "concrete_walk2.ogg", BaseStepPath + "concrete_walk3.ogg", BaseStepPath + "concrete_walk4.ogg",
										 BaseStepPath + "concrete_walk5.ogg", BaseStepPath + "concrete_walk6.ogg", BaseStepPath + "concrete_walk7.ogg", BaseStepPath + "concrete_walk8.ogg" };
			//
			SoundPathsStep["SqueakyWood"] = { BaseStepPath + "squeakywood_walk1.ogg", BaseStepPath + "squeakywood_walk2.ogg", BaseStepPath + "squeakywood_walk3.ogg", BaseStepPath + "squeakywood_walk4.ogg",
										 BaseStepPath + "squeakywood_walk5.ogg", BaseStepPath + "squeakywood_walk6.ogg", BaseStepPath + "squeakywood_walk7.ogg", BaseStepPath + "squeakywood_walk8.ogg" };


			// Break/Place (Modify)
			std::string BaseModifyPath = "Res/Sounds/Default/Modify/";
			SoundPathsModify["Cloth"] = { BaseModifyPath + "cloth1.ogg", BaseModifyPath + "cloth2.ogg", BaseModifyPath + "cloth3.ogg", BaseModifyPath + "cloth4.ogg" };

			SoundPathsModify["Grass"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };
			SoundPathsModify["Dirt"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };
			SoundPathsModify["Leaves"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };

			SoundPathsModify["Gravel"] = { BaseModifyPath + "gravel1.ogg", BaseModifyPath + "gravel2.ogg", BaseModifyPath + "gravel3.ogg", BaseModifyPath + "gravel4.ogg" };
			SoundPathsModify["Sand"] = { BaseModifyPath + "sand1.ogg", BaseModifyPath + "sand2.ogg", BaseModifyPath + "sand3.ogg", BaseModifyPath + "sand4.ogg" };
			SoundPathsModify["Snow"] = { BaseModifyPath + "snow1.ogg", BaseModifyPath + "snow2.ogg", BaseModifyPath + "snow3.ogg", BaseModifyPath + "snow4.ogg" };
			SoundPathsModify["Stone"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
			SoundPathsModify["Concrete"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
			SoundPathsModify["Metal"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
			SoundPathsModify["Wood"] = { BaseModifyPath + "wood1.ogg", BaseModifyPath + "wood2.ogg", BaseModifyPath + "wood3.ogg", BaseModifyPath + "wood4.ogg" };
			SoundPathsModify["SqueakyWood"] = { BaseModifyPath + "wood1.ogg", BaseModifyPath + "wood2.ogg", BaseModifyPath + "wood3.ogg", BaseModifyPath + "wood4.ogg" };
			SoundPathsModify["Marble"] = { BaseStepPath + "marble1.ogg", BaseStepPath + "marble2.ogg", BaseStepPath + "marble3.ogg", BaseStepPath + "marble4.ogg" };
		}

		else {
			std::string BaseStepPath = "Res/Sounds/Alt/Step/";
			SoundPathsStep["Cloth"] = { BaseStepPath + "cloth1.ogg", BaseStepPath + "cloth2.ogg", BaseStepPath + "cloth3.ogg", BaseStepPath + "cloth4.ogg" };
			SoundPathsStep["Grass"] = { BaseStepPath + "grass1.ogg", BaseStepPath + "grass2.ogg", BaseStepPath + "grass3.ogg", BaseStepPath + "grass5.ogg" };
			SoundPathsStep["Dirt"] = { BaseStepPath + "grass1.ogg", BaseStepPath + "grass2.ogg", BaseStepPath + "grass3.ogg", BaseStepPath + "grass5.ogg" };
			SoundPathsStep["Leaves"] = { BaseStepPath + "grass1.ogg", BaseStepPath + "grass2.ogg", BaseStepPath + "grass3.ogg", BaseStepPath + "grass5.ogg" };
			SoundPathsStep["Gravel"] = { BaseStepPath + "gravel1.ogg", BaseStepPath + "gravel2.ogg", BaseStepPath + "gravel3.ogg", BaseStepPath + "gravel4.ogg" };
			SoundPathsStep["Sand"] = { BaseStepPath + "sand1.ogg", BaseStepPath + "sand2.ogg", BaseStepPath + "sand3.ogg", BaseStepPath + "sand4.ogg" };
			SoundPathsStep["Snow"] = { BaseStepPath + "snow1.ogg", BaseStepPath + "snow2.ogg", BaseStepPath + "snow3.ogg", BaseStepPath + "snow4.ogg" };
			SoundPathsStep["Stone"] = { BaseStepPath + "stone1.ogg", BaseStepPath + "stone2.ogg", BaseStepPath + "stone3.ogg", BaseStepPath + "stone5.ogg" };
			SoundPathsStep["Concrete"] = { BaseStepPath + "stone1.ogg", BaseStepPath + "stone2.ogg", BaseStepPath + "stone3.ogg", BaseStepPath + "stone5.ogg" };
			SoundPathsStep["Wood"] = { BaseStepPath + "wood1.ogg", BaseStepPath + "wood2.ogg", BaseStepPath + "wood3.ogg", BaseStepPath + "wood4.ogg" };
			SoundPathsStep["SqueakyWood"] = { BaseStepPath + "wood1.ogg", BaseStepPath + "wood2.ogg", BaseStepPath + "wood3.ogg", BaseStepPath + "wood4.ogg" };
			SoundPathsStep["Metal"] = { BaseStepPath + "metal1.ogg", BaseStepPath + "metal2.ogg", BaseStepPath + "metal3.ogg", BaseStepPath + "metal4.ogg" };
			SoundPathsStep["Marble"] = { BaseStepPath + "marble1.ogg", BaseStepPath + "marble2.ogg", BaseStepPath + "marble3.ogg", BaseStepPath + "marble4.ogg" };

			// Break/Place (Modify)
			std::string BaseModifyPath = "Res/Sounds/Alt/Modify/";
			SoundPathsModify["Cloth"] = { BaseModifyPath + "cloth1.ogg", BaseModifyPath + "cloth2.ogg", BaseModifyPath + "cloth3.ogg", BaseModifyPath + "cloth4.ogg" };
			SoundPathsModify["Grass"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };
			SoundPathsModify["Dirt"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };
			SoundPathsModify["Leaves"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };
			SoundPathsModify["Gravel"] = { BaseModifyPath + "gravel1.ogg", BaseModifyPath + "gravel2.ogg", BaseModifyPath + "gravel3.ogg", BaseModifyPath + "gravel4.ogg" };
			SoundPathsModify["Sand"] = { BaseModifyPath + "sand1.ogg", BaseModifyPath + "sand2.ogg", BaseModifyPath + "sand3.ogg", BaseModifyPath + "sand4.ogg" };
			SoundPathsModify["Snow"] = { BaseModifyPath + "snow1.ogg", BaseModifyPath + "snow2.ogg", BaseModifyPath + "snow3.ogg", BaseModifyPath + "snow4.ogg" };
			SoundPathsModify["Stone"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
			SoundPathsModify["Concrete"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
			SoundPathsModify["Metal"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
			SoundPathsModify["Wood"] = { BaseModifyPath + "wood1.ogg", BaseModifyPath + "wood2.ogg", BaseModifyPath + "wood3.ogg", BaseModifyPath + "wood4.ogg" };
			SoundPathsModify["SqueakyWood"] = { BaseModifyPath + "wood1.ogg", BaseModifyPath + "wood2.ogg", BaseModifyPath + "wood3.ogg", BaseModifyPath + "wood4.ogg" };
			SoundPathsModify["Marble"] = { BaseStepPath + "marble1.ogg", BaseStepPath + "marble2.ogg", BaseStepPath + "marble3.ogg", BaseStepPath + "marble4.ogg" };

		}
	}
	
}