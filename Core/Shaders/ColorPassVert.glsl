#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec3 v_RayDirection;
out vec3 v_RayOrigin;
out vec2 v_TexCoords;
out vec2 v_CloudSourceOcclusion; // For basic specular occlusion


uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;


uniform sampler2D u_CloudProjectionVert;
uniform vec3 u_StrongerLightDirectionVert;
uniform vec3 u_SunDirectionVert;
uniform vec3 u_MoonDirectionVert;


vec2 ProjectDirection(vec3 Direction, vec2 TextureSize) 
{
    //vec2 TextureSize = vec2(textureSize(u_DebugTexture,0).xy);
	float TileSize = min(floor(TextureSize.x * 0.5f) / 1.5f, floor(TextureSize.y * 0.5f));
	float TileSizeDivided = (0.5f * TileSize) - 1.5f;
	vec2 CurrentCoordinate;

	if (abs(Direction.x) > abs(Direction.y) && abs(Direction.x) > abs(Direction.z)) 
    {
		Direction /= max(abs(Direction.x), 0.001f);
		CurrentCoordinate.x = Direction.y * TileSizeDivided + TileSize * 0.5f;
		CurrentCoordinate.y = Direction.z * TileSizeDivided + TileSize * (Direction.x < 0.0f ? 0.5f : 1.5f);
	} 
    
    else if (abs(Direction.y) > abs(Direction.x) && abs(Direction.y) > abs(Direction.z))
    {
		Direction /= max(abs(Direction.y), 0.001f);
		CurrentCoordinate.x = Direction.x * TileSizeDivided + TileSize * 1.5f;
		CurrentCoordinate.y = Direction.z * TileSizeDivided + TileSize * (Direction.y < 0.0f ? 0.5f : 1.5f);
	} 
    
    else 
    {
		Direction /= max(abs(Direction.z), 0.001f);
		CurrentCoordinate.x = Direction.x * TileSizeDivided + TileSize * 2.5f;
		CurrentCoordinate.y = Direction.y * TileSizeDivided + TileSize * (Direction.z < 0.0f ? 0.5f : 1.5f);
	}

	return CurrentCoordinate / max(TextureSize, 0.01f);
}

void main()
{
	const vec2 bayerSequenceOffsets[16] = vec2[16](vec2(0, 3) / 16.0, vec2(8, 11) / 16.0, vec2(2, 1) / 16.0, vec2(10, 9) / 16.0, vec2(12, 15) / 16.0, vec2(4, 7) / 16.0, vec2(14, 13) / 16.0, vec2(6, 5) / 16.0, vec2(3, 0) / 16.0, vec2(11, 8) / 16.0, vec2(1, 2) / 16.0, vec2(9, 10) / 16.0, vec2(15, 12) / 16.0, vec2(7, 4) / 16.0, vec2(13, 14) / 16.0, vec2(5, 6) / 16.0);

	gl_Position = vec4(a_Position.xy, 1.0f, 1.0f);
	v_TexCoords.xy = a_TexCoords.xy;
	vec2 Position = a_Position;
	vec4 clip = vec4(Position.xy, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	v_RayDirection = vec3(u_InverseView * eye);
	v_RayOrigin = u_InverseView[3].xyz;
	
	//vec2 Dim = textureSize(u_CloudProjectionVert, 0).xy;
	//v_CloudSourceOcclusion = vec2(texture(u_CloudProjectionVert, ProjectDirection(u_SunDirectionVert, Dim)).w, texture(u_CloudProjectionVert, ProjectDirection(u_SunDirectionVert, Dim)).w);
	v_CloudSourceOcclusion = vec2(0.0001f);								
}