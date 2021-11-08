/* VoxelRT - A Voxel Raytracing Engine
Written by : Samuel Rasquinha

Contributors (In no particular order.) : 
Lars
Snurf
Fuzdex (ShadaxStack)
Moonsheep
Tui Vao (Tester)
UglySwedishFish
Telo

Resources : 
https://github.com/BrutPitt/glslSmartDeNoise/
https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
ScratchAPixel
https://google.github.io/filament/Filament.md
http://advances.realtimerendering.com/s2015/The%20Real-time%20Volumetric%20Cloudscapes%20of%20Horizon%20-%20Zero%20Dawn%20-%20ARTR.pdf
http://www.diva-portal.org/smash/get/diva2:1223894/FULLTEXT01.pdf

Textures : 
Most are from the RealismMats texture pack and umsoea. The rest are from CC0
I do not make any money from this project. This project is purely educational

Sounds : 
- Taken from the minecraft soundpack "Enhanced sounds" 
- 2 Other sounds taken from websites, CC0 licensed

-- Please do not copy my work and claim it as your own --
*/ 

///////////////////////////////////////////////////////////////////////////////////////////
//    																					 //
//    MIT License																		 //
//    																					 //
//    Copyright (c) 2021 Samuel Rasquinha												 //
//    																					 //
//    Permission is hereby granted, free of charge, to any person obtaining a copy		 //
//    of this software and associated documentation files (the "Software"), to deal		 //
//    in the Software without restriction, including without limitation the rights		 //
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell			 //
//    copies of the Software, and to permit persons to whom the Software is				 //
//    furnished to do so, subject to the following conditions:							 //
//    																					 //
//    The above copyright notice and this permission notice shall be included in all	 //
//    copies or substantial portions of the Software.									 //
//    																					 //
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR		 //
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,			 //
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE		 //
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER			 //
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,		 //
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE		 //
//    SOFTWARE.																			 //
//    																		             //
///////////////////////////////////////////////////////////////////////////////////////////

/* Controls 
F1 -> Lock/Unlock mouse
W, S, A, D -> Move 
Space -> Jump/Fly
Shift -> Accelerate down
Q/E -> Change current block
F -> Toggle freefly
C -> Toggle Collisions (only toggles IF on freefly mode) 
ESC -> Save and quit
V -> Toggle VSync
F2 -> Recompile shaders
ImGui Windows -> 
	Window 1 : Various Other Settings (Resolution settings) 
	Window 2 : Player sensitivity, speed and sound options 
*/

/*
Thankyou.
*/

// Main Pipeline
#include "Core/Pipeline.h"

// Main
int main()
{
	VoxelRT::MainPipeline::StartPipeline();
	return 0;
}

//
