if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Project Setup
cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
project(VoxelRT)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-fno-rtti -pedantic-errors -fomit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -static")

# GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(ext/glfw-3.3.2)

set(CMAKE_CXX_FLAGS "-fno-rtti -pedantic-errors -fomit-frame-pointer -mavx2 -mavx512f -msse2 -msse4.1")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -static")

add_executable(VoxelRT
        Core/AORenderer.cpp
        Core/AORenderer.h
		Core/AtmosphereRenderCubemap.cpp
        Core/AtmosphereRenderCubemap.h
		Core/AtmosphereRenderer.cpp
        Core/AtmosphereRenderer.h
        Core/Audio.h
        Core/Block.h
        Core/BlockDatabase.h
        Core/BlockDatabase.cpp
		Core/BlockDatabaseParser.h
        Core/BlockDatabaseParser.cpp
		Core/BlockDataSSBO.h
        Core/BlockDataSSBO.cpp
		Core/BloomFBO.h
        Core/BloomFBO.cpp
		Core/BloomRenderer.h
        Core/BloomRenderer.cpp
		Core/BlueNoiseDataSSBO.h
        Core/BlueNoiseDataSSBO.cpp
		Core/ColorPassFBO.h
        Core/ColorPassFBO.cpp
		Core/FpsCamera.h
        Core/FpsCamera.cpp
		Core/Macros.h
		Core/OrthographicCamera.h
        Core/OrthographicCamera.cpp
		Core/Particle.h
		Core/ParticleSystem.h
        Core/ParticleSystem.cpp
		Core/Pipeline.h
        Core/Pipeline.cpp
		Core/Player.h
        Core/Player.cpp
		Core/Renderer2D.h
        Core/Renderer2D.cpp
		Core/RenderTexture.h
        Core/RenderTexture.cpp
		Core/ShaderManager.h
        Core/ShaderManager.cpp
		Core/SoundManager.h
        Core/SoundManager.cpp
		Core/Texture3D.h
        Core/Texture3D.cpp
		Core/World.h
        Core/World.cpp
		Core/WorldFileHandler.h
        Core/WorldFileHandler.cpp
		Core/WorldGenerator.h
        Core/WorldGenerator.cpp
		
		Core/Application/Application.h
		Core/Application/Application.cpp
		Core/Application/Logger.h
		Core/Application/Logger.cpp
		
		Core/Clouds/CloudFBO.h
		Core/Clouds/CloudFBO.cpp
		Core/Clouds/CloudRenderer.h
		Core/Clouds/CloudRenderer.cpp
		Core/Clouds/NoiseRenderer.h
		Core/Clouds/NoiseRenderer.cpp
		Core/Clouds/NoiseTexture3D.h
		Core/Clouds/NoiseTexture3D.cpp
		
		Core/GLClasses/ComputeShader.h
		Core/GLClasses/ComputeShader.cpp
		Core/GLClasses/CubeTextureMap.h
		Core/GLClasses/CubeTextureMap.cpp
		Core/GLClasses/Fps.h
		Core/GLClasses/Fps.cpp
		Core/GLClasses/Framebuffer.h
		Core/GLClasses/Framebuffer.cpp
		Core/GLClasses/FramebufferRed.h
		Core/GLClasses/FramebufferRed.cpp
		Core/GLClasses/IndexBuffer.h
		Core/GLClasses/IndexBuffer.cpp
		Core/GLClasses/stb_image.h
		Core/GLClasses/stb_image.cpp
		Core/GLClasses/Texture.h
		Core/GLClasses/Texture.cpp
		Core/GLClasses/TextureArray.h
		Core/GLClasses/TextureArray.cpp
		Core/GLClasses/VertexArray.h
		Core/GLClasses/VertexArray.cpp
		Core/GLClasses/VertexBuffer.h
		Core/GLClasses/VertexBuffer.cpp

		Core/Physics/AABB.h
		
		Core/GLClasses/Random.h
		Core/GLClasses/Timer.h
       )

target_link_libraries(VoxelRT 
	glfw
	glad
	glm
	imgui::imgui
	irrklang)