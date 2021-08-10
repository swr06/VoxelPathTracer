# Voxel Path Tracer
A Voxel Path Tracing Engine implemented using C++ and the modern opengl programmable pipeline

# Features Implemented
- Voxel Ray Casting 
- Manhattan Distance Field acceleration (Distance fields are generated using compute shaders) 
- Procedural Terrain Generation (Plains, superflat) 
- Indirect Light Tracing (Diffuse, Specular, AO)
- Temporal Indirect Light Filtering (Specular reprojection done using hit distance.) 
- Parallax Occlusion Mapping (Buggy with materials adapted from minecraft texture packs!)
- Normal mapping
- PBR Texture support 
- Physically based lighting (Cook torrance BRDF)
- Ray traced contact hardening soft shadows (with shadow denoising)
- Flexible material system
- Temporal Anti Aliasing
- FXAA
- Tonemapping, Gamma correction
- Ray traced reflections 
- Ray traced rough reflections (Importance Sampled GGX)
- Naive world saving and loading
- Physically based atmosphere 
- God Rays (Screen space) 
- Per Pixel Emissive Materials
- Spatial filtering (Atrous and SVGF)
- Screen Space Ambient Occlusion
- Ray Traced Ambient Occlusion
- Bloom (Mip based) 
- Volumetric clouds
- Particle system
- Spherical harmonics (to encode indirect radiance data, used for both specular and diffuse) 
- Alpha testing (Slightly buggy with the diffuse indirect :/) 
- Fully 3D Audio

# Features to implement
- Foliage
- Refractions (Will be done in screen space)
- Glass rendering 
- Water Rendering (Tesselation with FFT? tessendorf waves?)

# Performance

- 30 FPS on a Vega 8 iGPU on the default settings
- 180 - 200 FPS on a GTX 1080Ti

# Note
- The path tracer is still HEAVILY WIP, the current state is not a representation of the final version.
- The path tracer has been tested on AMD Vega iGPUs, AMD GPUs, Nvidia pascal, turing and ampere cards.
- The path tracer is *not* guarenteed to work on ANY Intel GPUs
- It needs OpenGL 4.5 (Uses compute shaders and other features from OpenGL 4.5), if the window fails to initialize, then your GPU does not support the required OpenGL version (OpenGL 4.5) 

# Resources used
- https://github.com/BrutPitt/glslSmartDeNoise/
- https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
- ScratchAPixel
- https://google.github.io/filament/Filament.md
- http://advances.realtimerendering.com/s2015/The%20Real-time%20Volumetric%20Cloudscapes%20of%20Horizon%20-%20Zero%20Dawn%20-%20ARTR.pdf
- http://www.diva-portal.org/smash/get/diva2:1223894/FULLTEXT01.pdf
- https://github.com/NVIDIA/Q2RTX
- https://cg.ivd.kit.edu/publications/2017/svgf/svgf_preprint.pdf
- https://teamwisp.github.io/research/svfg.html
- Textures from "Realism Mats", "Umsoea", CC0 textures and Quixel Megascans 

# Thanks
- [Fuzdex](https://github.com/Shadax-stack)
- [UglySwedishFish](https://github.com/UglySwedishFish)
- [Lars](https://github.com/Ciwiel3/)
- [Snurf](https://github.com/AntonioFerreras)
- [Telo](https://github.com/StormCreeper)
- [Tui Vao](https://github.com/Tui-Vao)
- [Moonsheep](https://github.com/jlagarespo)

# License
- MIT

# Notice
This project is purely educational. I own none of the assets. All the rights go to their respective owners.
