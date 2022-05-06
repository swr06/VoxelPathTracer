# Voxel Path Tracer
An ***EXPERIMENTAL*** and ***UNFINISHED*** Voxel Path Tracing Engine which has an emphasis on performance and graphics. This engine was mostly made as an experiment and a tool to help me learn more about light transport, physically based rendering, volumetrics and intersection algorithms.
This engine was coded from scratch using C++17 and the modern OpenGL programmable pipeline (OpenGL 4.5+ and GLSL version 430+ required).

## Project status 
This project has (mostly) been abandoned and new *big* features will not be worked on anytime soon.
I intend to rewrite a voxel tracing engine in the future with much cleaner code/syntax and with newer and better acceleration structures so that larger and more detailed worlds can be rendered performantly!

## Features Implemented

### Rendering 
- Voxel Ray Casting 
- Manhattan distance field acceleration structure (Distance field generated using compute shaders on the GPU) 
- Deferred rendering pipeline

### Lighting 
- Direct lighting based on the Cook torrance BRDF
- Path traced lighting (Direct shadows, Global Illumination, Rough/Smooth reflections and Ambient Occlusion)
- Spherical harmonic projection for indirect lighting 
- (Approximate) Screenspace Subsurface Scattering
- Dynamic atmosphere/sky rendering
- Relief Parallax Mapping

### Denoising
- Temporal Denoiser/Reprojection (used for Pathtraced Lighting, Volumetric Clouds and Antialiasing)
- Screenspace Spatial Denoiser (SVGF, Atrous, Gaussian and other specialized denoisers for direct shadow/reflections)

### Post Process

- AO : SSAO
- Misc : DOF, Bloom, Lens flare, Chromatic aberration, Basic SS god rays, night stars and nebula
- Image : CAS
- Anti Aliasing : FXAA, TXAA/TAA-U
- Other Effects : Tonemapping, Gamma correction, color dithering, color grading, purkinje effect, film grain

### Volumetrics 
- 3D Volumetric Clouds (+ 2D cirrus cloud layer)
- LPV for volumetrics from light sources (could be used for distant lighting or performant second bounce gi in the future)

### Others
- Basic Particle system
- 3D Sound
- Material/Block database system (with around 100 block types with 512x512 textures (Albedo, normals, metalness, roughness, indirect occlusion & displacement))
- Procedural Terrain Generation (Plains, superflat) 
- World saving/loading
- Minecraft world loading 

## Half Implemented / WIP Features
- Histogram based auto exposure
- Anisotropic raytraced reflections (for materials such as brushed metal, based on Disney Principled BRDF)

## Todo Features / QOL Improvements
- Glass/water with stochastic refractions
- Worldspace MIS 

## Performance Metrics 

- 20 - 22 FPS on a Vega 8 iGPU on the default low preset. (@ 768p)
- 30 - 34 FPS on a Vega 11 iGPU on the default low preset. (@ 768p)
- 150 - 160 FPS on a GTX 1080Ti. (@ 1080p)

## Note
- This project is still very experimental, it is not a finished product as of yet.
- This engine has *NOT* been tested on Intel GPUs/iGPUs.
- This engine requires OpenGL 4.5, if the window fails to initialize, then there is a good chance that your GPU does not support the required OpenGL version.
- If you want to report an issue/bug, then you can contact me on discord or, alternatively, via email. (See github profile page)
- See `Controls.txt` for the controls (Or look at the console when you start up the program.)

## Credits (Testing, programming and understanding)
- [Lars](https://github.com/Ciwiel3/)
- [UglySwedishFish](https://github.com/UglySwedishFish)
- [Snurf](https://github.com/AntonioFerreras)
- [Fuzdex](https://github.com/Shadax-stack)
- [Telo](https://github.com/StormCreeper)
- [Moonsheep](https://github.com/jlagarespo)
- [Tui Vao](https://github.com/Tui-Vao)

## License
- MIT (See LICENSE)

## Notice
This project is purely educational. I own none of the assets. All the rights go to their respective owners.

## Screenshots 

</br>

![glare](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/glare.png)

</br>

</br>

![day](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/day.png)

</br>

![day2](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/day2.png)

</br>

![day3](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/day3.png)

</br>

![dof](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/dof.png)

</br>

</br>

![terrain](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/terrain.png)

</br>

![night](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/night1.png)

</br>

![night2](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/night2.png)

</br>

![amogsus](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/amogus.png)

</br>

![amogsus2](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/amogus2.png)

</br>

![vol1](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/vol1.png)

</br>

![vol2](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/vol2.png)

</br>

![clood2](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/volclouds2.png)

</br>

![postprocesshahahahahhahahahALSOifyouSeeThisIWillGiveYouACookieAddMeOnDiscord](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/postprocessgobrr.png)

</br>

![shadowwww](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/shadow.png)

</br>

</br>

# Supporting

If you'd like to support this project, consider starring it. Thanks!
