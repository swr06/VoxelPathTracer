# Voxel Path Tracer
An ***UNFINISHED*** and ***VERY EXPERIMENTAL*** Voxel Path Tracing Engine which has an emphasis on performance and graphics. This engine was mostly made as an experiment and a tool to help me learn more about light transport, volumetrics and intersection algorithms.
This engine was coded from scratch using C++17 and the modern OpenGL programmable pipeline (OpenGL 4.5+ and GLSL version 430+ required)

## Project status 
This project has been abandoned and will not be worked on anytime soon (The projects that I work on regularly are usually kept private). 
I intend to rewrite a voxel tracing engine in the future with much cleaner code/syntax and with newer and better acceleration structures so that larger scale worlds can be rendered with ease!

## Features Implemented

### Rendering 
- Voxel Ray Casting 
- Manhattan Distance Field acceleration (Distance fields are generated using compute shaders) 
- Procedural Terrain Generation (Plains, superflat) 

### Lighting 
- Direct lighting based on the Cook torrance BRDF
- Path traced lighting (Direct shadows, Global Illumination, Rough/Smooth reflections, Ambient Occlusion)
- Screen space SSS
- Dynamic atmosphere rendering 

### Post Process

- AO : SSAO
- Misc : DOF, Bloom, Lens flare, Chromatic aberration, Basic SS god rays, night stars and nebula
- Image : CAS
- Anti Aliasing : FXAA, TXAA
- Color Compositing : Tonemapping, Gamma correction, color dithering, color grading, purkinje effect, film grain

### Volumetrics 
- 3D Volumetric Clouds
- LPV for volumetrics from light sources (could be used for distant lighting or quick second bounce gi in the future)

### Others
- Basic Particle system
- 3D Sound
- Material/Block database system (with around 100 block types with 512x512 textures (Albedo, normals, metalness, roughness, indirect occlusion & displacement))
- World saving/loading
- Minecraft world loading 

## Half Implemented / WIP Features
- Histogram based auto exposure
- Parallax mapping (Parallax Relief Mapping, Parallax Occlusion Mapping and Ground Truth Parallax)
- Anisotropic raytraced reflections (for materials such as brushed metal, based on principled disney BRDF)

## Todo Features / QOL Improvements
- Bug fixes
- Better Player Controller
- Foliage
- (Preferably sthochastic) Refractions
- Glass rendering 
- Water Rendering (Tesselation with FFT? tessendorf waves?)
- Raytraced realistic sound
- Weather?

## Performance Metrics 

- 24 FPS on a Vega 8 iGPU on the default settings.
- 180 - 200 FPS on a GTX 1080Ti.

## Note
- This project is still not finished, the current state is not a representation of the final version.
- It has been tested on AMD Vega iGPUs, AMD GPUs, Nvidia pascal, turing and ampere cards.
- It is *not* guarenteed to work on ANY Intel GPUs
- It needs OpenGL 4.5 (Uses compute shaders and other features from OpenGL 4.5), if the window fails to initialize, then your GPU does not support the required OpenGL version.
- If you want to report an issue/bug, then you can contact me on discord or, alternatively, via email. (See github profile)
- See `Controls.txt` for the controls (Or look at the console when you start up the program!)

## Credits (In no particular order.)
- [Fuzdex](https://github.com/Shadax-stack)
- [UglySwedishFish](https://github.com/UglySwedishFish)
- [Lars](https://github.com/Ciwiel3/)
- [Snurf](https://github.com/AntonioFerreras)
- [Telo](https://github.com/StormCreeper)
- [Tui Vao](https://github.com/Tui-Vao)
- [Moonsheep](https://github.com/jlagarespo)

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

![clood](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/volclouds.png)

</br>

![clood2](https://github.com/swr06/VoxelPathTracer/blob/Project-Main/Screenshots/volclouds2.png)

</br>





