#pragma once

#include <iostream>
#include <string>
#include <glad/glad.h>

#include "../GLClasses/Framebuffer.h"
#include "NoiseTexture3D.h"
#include "../GLClasses/Shader.h"
#include "../GLClasses/VertexArray.h"
#include "../GLClasses/VertexBuffer.h"
#include "../GLClasses/IndexBuffer.h"

namespace Clouds
{
	void RenderNoise(NoiseTexture3D& tex, int slices, bool detail = false);
	void RenderCurlNoise(GLClasses::Framebuffer& fbo);
}