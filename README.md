# Voxel Path Tracer
A Voxel Path Tracing Engine implemented using C++ and the modern opengl programmable pipeline

# Features Implemented
- Voxel Ray Casting
- Manhattan Distance Field acceleration (Distance fields are generated using compute shaders) 
- Procedural Terrain Generation
- Indirect Light Tracing 
- Temporal Indirect Light Filtering
- Parallax Occlusion Mapping
- Normal mapping
- PBR Texture support
- Physically based lighting
- Ray traced contact hardening soft shadows
- Flexible material system
- Temporal Anti Aliasing
- FXAA
- Tonemapping, Gamma correction
- Ray traced reflections 
- Ray traced rough reflections (Importance Sampled GGX)
- Checkerboard raytracing 
- Naive world saving and loading
- Physically based atmosphere 
- God Rays (Screen space) 
- Emissive Materials
- Atrous Spatial filtering
- Screen Space Ambient Occlusion
- Ray Traced Ambient Occlusion
- Lens Flare
- Bloom (Mip based) 
- Volumetric clouds
- Particle system

# Features to implement
- Spherical Harmonics to increase normal map sharpness
- Alpha testing
- Foliage
- Refractions (Will be done in screen space)
- Glass rendering 
- Water Rendering (Tesselation with FFT? tessendorf waves?)

# Performance

- 30 FPS on a Vega 8 iGPU on the default settings
- 180 - 200 FPS on a GTX 1080Ti

# Note
- The path tracer has been tested on AMD Vega iGPUs, AMD GPUs, Nvidia pascal, turing and ampere cards.
- It needs OpenGL 4.5 (Uses compute shaders and other features from OpenGL 4.5), if the window fails to initialize, then your GPU does not support the required OpenGL version (OpenGL 4.5) 

# Resources used
- https://github.com/BrutPitt/glslSmartDeNoise/
- https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
- ScratchAPixel
- https://google.github.io/filament/Filament.md
- http://advances.realtimerendering.com/s2015/The%20Real-time%20Volumetric%20Cloudscapes%20of%20Horizon%20-%20Zero%20Dawn%20-%20ARTR.pdf
- http://www.diva-portal.org/smash/get/diva2:1223894/FULLTEXT01.pdf

# Thanks
- [Fuzdex](https://github.com/Shadax-stack)
- [UglySwedishFish](https://github.com/UglySwedishFish)
- [Lars](https://github.com/Ciwiel3/)
- Snurf (Founder of ApolloRT)
- [Telo](https://github.com/StormCreeper)
- [Tui Vao](https://github.com/Tui-Vao)
- [Moonsheep](https://github.com/jlagarespo)

# License
- MIT

# Notice
This project is purely educational. I own none of the assets. All the rights go to their respective owners.
