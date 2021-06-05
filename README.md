# Voxel Tracing
A Voxel Ray Tracing Engine implemented using C++ and the modern opengl programmable pipeline

# Features Implemented
- Voxel Ray Casting
- Procedural Terrain Generation
- Diffuse Tracing 
- Temporal Diffuse Filtering
- Normal mapping
- PBR Texture support
- Physically based lighting
- Fast Shadow Tracing
- Flexible material system
- Temporal Anti Aliasing
- FXAA
- Tonemapping, Gamma correction
- Bilateral Upsampling (For SSAO, Trace)
- Ray traced reflections 
- Ray traced rough reflections (Importance Sampled GGX)
- Naive world saving and loading
- Alpha testing (Shadows and tracing)
- Atmosphere 
- God Rays (Screen space) 
- Emissive Materials
- Basic spatial filtering (Not depth or normal aware as of yet.)
- Screen Space Ambient Occlusion
- Ray Traced Ambient Occlusion
- Lens Flare
- Bloom (Mip based) 

# Features to implement
- Checkerboarding for the raw diffuse trace
- Better Soft Shadows
- Refractions (Will be done in screen space)
- Glass rendering 
- Better denoising
- Water Rendering (Tesselation with FFT? tessendorf waves?)
- Voxel intersection acceleration structure (DAGs? Octrees?) 
- Parallax Occlusion Mapping (Try tracing through the volume, Maybe?)
- Cheap clouds 

# Performance

- 30 FPS on a Vega 8 iGPU on the default settings

# Resources used
- https://github.com/BrutPitt/glslSmartDeNoise/
- https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
- ScratchAPixel
- https://google.github.io/filament/Filament.md

# Thanks
- Fuzdex (ShadaxStack)
- UglySwedishFish
- Snurf (Founder of ApolloRT)
- Telo 
- Moonsheep

# Notice
This project is purely educational. I own none of the assets. All the rights go to their respective owners.
