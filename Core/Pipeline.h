#pragma once


#include <iostream>

#include "Application/Application.h"
#include "GLClasses/VertexArray.h"
#include "GLClasses/VertexBuffer.h"
#include "GLClasses/IndexBuffer.h"
#include "GLClasses/Shader.h"
#include "FpsCamera.h"
#include "GLClasses/Fps.h"
#include "World.h"
#include "GLClasses/Framebuffer.h"
#include "WorldGenerator.h"
#include "GLClasses/CubeTextureMap.h"
#include "BlockDatabase.h"
#include "Player.h"
#include "GLClasses/FramebufferRed.h"
#include "ColorPassFBO.h"
#include "GLClasses/Texture.h"
#include "Renderer2D.h"
#include "WorldFileHandler.h"
#include "AtmosphereRenderer.h"
#include "AtmosphereRenderCubemap.h"
#include "BloomFBO.h"
#include "BloomRenderer.h"
#include "Clouds/CloudRenderer.h"

namespace VoxelRT
{
	namespace MainPipeline
	{
		void StartPipeline();
	}
}