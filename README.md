# Voxel Path Tracer
A Voxel Path Tracing Engine which has an emphasis on performance and graphics. 
This engine was implemented using C++ and the modern opengl programmable pipeline.

# Features Implemented
- Voxel Ray Casting 
- Over 100+ block types with 512x512 textures (Albedo, normals, metalness, roughness, indirect occlusion & displacement)
- Manhattan Distance Field acceleration (Distance fields are generated using compute shaders) 
- Procedural Terrain Generation (Plains, superflat) 
- Indirect Light Tracing (Diffuse, Specular, AO)
- Temporal Indirect Light Filtering (Specular reprojection done using hit distance.) 
- Parallax Occlusion Mapping (Buggy with materials adapted from minecraft texture packs!)
- Normal mapping
- PBR Texture support 
- Physically based lighting (Cook torrance BRDF)
- Ray traced contact hardening soft shadows (with specialized shadow denoiser)
- Flexible material system
- Temporal Anti Aliasing (And temporal supersampling)
- FXAA
- HDR Tonemapping, Gamma correction
- Post Processing (Lens flare, Purkinje Effect, Film grain, Color Dither, Color grading with an LUT, Chromatic Aberration, Screen Space God Rays etc.)
- Ray traced reflections 
- Ray traced rough reflections (Importance Sampled GGX with specialized specular denoiser)
- Naive world saving and loading
- Physically based atmosphere
- Per Pixel Emissive Materials
- Spatial filtering (Atrous and SVGF)
- Spatial Upscaling
- Screen Space Ambient Occlusion 
- Ray Traced Ambient Occlusion
- Bloom (Mip based) 
- Volumetric 3D clouds (with temporal and spatial upscaling/filtering)
- Particle system
- Spherical harmonics (Second order SH, used to encode indirect radiance data, used for both specular and diffuse) 
- Alpha testing
- Fully 3D Audio
- Light Propogation Volumes with custom smooth interpolation and support for colored lights
- Point Light Volumetrics using ray marching and the Light Propogation Volume
- Contrast Adaptive Sharpening
- Auto exposure based on a luminance histogram (WIP.)

# Features to implement (Or atleast.. in mind)
- Foliage
- Sunsets 
- (Preferably sthochastic) Refractions
- Glass rendering 
- Water Rendering (Tesselation with FFT? tessendorf waves?)
- Raytraced realistic sound

# Performance

- 30 FPS on a Vega 8 iGPU on the default settings
- 180 - 200 FPS on a GTX 1080Ti

# Note
- This project is still not finished, the current state is not a representation of the final version.
- It has been tested on AMD Vega iGPUs, AMD GPUs, Nvidia pascal, turing and ampere cards.
- It is *not* guarenteed to work on ANY Intel GPUs
- It needs OpenGL 4.5 (Uses compute shaders and other features from OpenGL 4.5), if the window fails to initialize, then your GPU does not support the required OpenGL version.
- I don't work on this project a lot anymore, newer features might be delayed.
- If you want to report an issue, then you can contact me on discord (swr#1793)
- See `Controls.txt` for the controls (Or look at the console when you start up the program!) 
- Amogus.

# Resource List
- https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf (DDA)
- ScratchAPixel
- https://google.github.io/filament/Filament.md
- http://advances.realtimerendering.com/s2015/The%20Real-time%20Volumetric%20Cloudscapes%20of%20Horizon%20-%20Zero%20Dawn%20-%20ARTR.pdf
- http://www.diva-portal.org/smash/get/diva2:1223894/FULLTEXT01.pdf
- https://github.com/NVIDIA/Q2RTX
- https://cg.ivd.kit.edu/publications/2017/svgf/svgf_preprint.pdf
- https://www.youtube.com/watch?v=2XXS5UyNjjU (A Talk on TAA, from Inside.)
- https://teamwisp.github.io/research/svfg.html
- http://magnuswrenninge.com/wp-content/uploads/2010/03/Wrenninge-OzTheGreatAndVolumetric.pdf
- https://media.contentapi.ea.com/content/dam/eacom/frostbite/files/s2016-pbs-frostbite-sky-clouds-new.pdf
- Textures from CC0 textures (or AmbientCG), Quixel Megascans and textures.com
- Exactly 5 block textures are taken from the realism mats texture pack and 5 more from the umsoea texture pack. And 3 from the patrix texture pack. (If you want them removed, I will happily do so)

# Thanks (In no particular order.)
- [Fuzdex](https://github.com/Shadax-stack)
- [UglySwedishFish](https://github.com/UglySwedishFish)
- [Lars](https://github.com/Ciwiel3/)
- [Snurf](https://github.com/AntonioFerreras)
- [Telo](https://github.com/StormCreeper)
- [Tui Vao](https://github.com/Tui-Vao)
- [Moonsheep](https://github.com/jlagarespo)

# License
- MIT (See LICENSE)

# Notice
This project is purely educational. I own none of the assets. All the rights go to their respective owners.
