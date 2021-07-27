#include "SoundManager.h"

#include "Audio.h"

extern float VoxelRT_VolumeMultiplier;

namespace VoxelRT
{
	static irrklang::ISoundEngine* MainSoundEngine;
	static std::unordered_map<std::string, std::array<std::string, 4>> SoundPathsStep;
	static std::unordered_map<std::string, std::array<std::string, 4>> SoundPathsModify;

	void SoundManager::InitializeSoundManager()
	{
		MainSoundEngine = irrklang::createIrrKlangDevice();

		if (!MainSoundEngine)
		{
			throw "Sound engine couldnt be initialized :/";
		}

		// Step :
		std::string BaseStepPath = "Res/Sounds/Step/";
		SoundPathsStep["Cloth"] = { BaseStepPath + "cloth1.ogg", BaseStepPath + "cloth2.ogg", BaseStepPath + "cloth3.ogg", BaseStepPath + "cloth4.ogg" };
		SoundPathsStep["Grass"] = { BaseStepPath + "grass1.ogg", BaseStepPath + "grass2.ogg", BaseStepPath + "grass3.ogg", BaseStepPath + "grass5.ogg" };
		SoundPathsStep["Gravel"] = { BaseStepPath + "gravel1.ogg", BaseStepPath + "gravel2.ogg", BaseStepPath + "gravel3.ogg", BaseStepPath + "gravel4.ogg" };
		SoundPathsStep["Sand"] = { BaseStepPath + "sand1.ogg", BaseStepPath + "sand2.ogg", BaseStepPath + "sand3.ogg", BaseStepPath + "sand4.ogg" };
		SoundPathsStep["Snow"] = { BaseStepPath + "snow1.ogg", BaseStepPath + "snow2.ogg", BaseStepPath + "snow3.ogg", BaseStepPath + "snow4.ogg" };
		SoundPathsStep["Stone"] = { BaseStepPath + "stone1.ogg", BaseStepPath + "stone2.ogg", BaseStepPath + "stone3.ogg", BaseStepPath + "stone5.ogg" };
		SoundPathsStep["Wood"] = { BaseStepPath + "wood1.ogg", BaseStepPath + "wood2.ogg", BaseStepPath + "wood3.ogg", BaseStepPath + "wood4.ogg" };
		SoundPathsStep["Metal"] = { BaseStepPath + "metal1.ogg", BaseStepPath + "metal2.ogg", BaseStepPath + "metal3.ogg", BaseStepPath + "metal4.ogg" };
		SoundPathsStep["Marble"] = { BaseStepPath + "marble1.ogg", BaseStepPath + "marble2.ogg", BaseStepPath + "marble3.ogg", BaseStepPath + "marble4.ogg" };
	
		// Break/Place (Modify)
		std::string BaseModifyPath = "Res/Sounds/Modify/";
		SoundPathsModify["Cloth"] = { BaseModifyPath + "cloth1.ogg", BaseModifyPath + "cloth2.ogg", BaseModifyPath + "cloth3.ogg", BaseModifyPath + "cloth4.ogg" };
		SoundPathsModify["Grass"] = { BaseModifyPath + "grass1.ogg", BaseModifyPath + "grass2.ogg", BaseModifyPath + "grass3.ogg", BaseModifyPath + "grass4.ogg" };
		SoundPathsModify["Gravel"] = { BaseModifyPath + "gravel1.ogg", BaseModifyPath + "gravel2.ogg", BaseModifyPath + "gravel3.ogg", BaseModifyPath + "gravel4.ogg" };
		SoundPathsModify["Sand"] = { BaseModifyPath + "sand1.ogg", BaseModifyPath + "sand2.ogg", BaseModifyPath + "sand3.ogg", BaseModifyPath + "sand4.ogg" };
		SoundPathsModify["Snow"] = { BaseModifyPath + "snow1.ogg", BaseModifyPath + "snow2.ogg", BaseModifyPath + "snow3.ogg", BaseModifyPath + "snow4.ogg" };
		SoundPathsModify["Stone"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
		SoundPathsModify["Metal"] = { BaseModifyPath + "stone1.ogg", BaseModifyPath + "stone2.ogg", BaseModifyPath + "stone3.ogg", BaseModifyPath + "stone4.ogg" };
		SoundPathsModify["Wood"] = { BaseModifyPath + "wood1.ogg", BaseModifyPath + "wood2.ogg", BaseModifyPath + "wood3.ogg", BaseModifyPath + "wood4.ogg" };
		SoundPathsModify["Marble"] = { BaseStepPath + "marble1.ogg", BaseStepPath + "marble2.ogg", BaseStepPath + "marble3.ogg", BaseStepPath + "marble4.ogg" };

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
		if (block <= 0)
		{
			return;
		}

		// 1 : step, 0 : modify
		std::string Type = type ? BlockDatabase::GetStepSound(block) : BlockDatabase::GetModifySound(block);
		auto& list = type ? SoundPathsStep : SoundPathsModify;

		if (list.find(Type) == list.end()) {
			std::cout << "\nCant play sound for block : " << block << "\n" ;
			return;
		}

		int rand_idx = rand() % 4;
		rand_idx = glm::clamp((int)rand_idx, (int)0, (int)3);
		std::string Path = list[Type][rand_idx];

		float v = type ? 6.0f : 6.0f; v *= VoxelRT_VolumeMultiplier;
		PlaySound(Path, p, type ? 4.0f : 8.0f, v, false);
	}

	void SoundManager::Destroy()
	{
		if (!MainSoundEngine) 
		{ 
			return; 
		}

		MainSoundEngine->drop();
	}
	
}