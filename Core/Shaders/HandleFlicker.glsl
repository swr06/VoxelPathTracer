// Basic anti flickering shader test
// doesnt work without reprojection
// at that point it becomes a temporal filter for the main image which I already have (TAA)

#version 330 core

layout (location = 0) out vec3 o_Result;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentTexture;
uniform sampler2D u_History[3];

const float FlickerThreshold = 0.35f;

float Luminance(vec3 x)
{
	return dot(x, vec3(0.2125f, 0.7154f, 0.0721f));
}

void main() {
	vec3 Current = texture(u_CurrentTexture, v_TexCoords).xyz;
	vec3 History_1 = texture(u_History[0], v_TexCoords).xyz;
	vec3 History_2 = texture(u_History[1], v_TexCoords).xyz;
	vec3 History_3 = texture(u_History[2], v_TexCoords).xyz;

	// gather Luminance
	float CurrentLuma = Luminance(Current);
	float HistoryLuma_1 = Luminance(History_1);
	float HistoryLuma_2 = Luminance(History_2);
	float HistoryLuma_3 = Luminance(History_3);

	// Use history to check for potential flickering :
	bool FlickerThere = (abs(CurrentLuma - HistoryLuma_1) > FlickerThreshold &&
						 abs(CurrentLuma - HistoryLuma_2) < FlickerThreshold) &&
						(abs(HistoryLuma_1 - HistoryLuma_2) > FlickerThreshold &&
						 abs(HistoryLuma_1 - HistoryLuma_3) < FlickerThreshold);

	// Mix
	vec3 Interpolated = Current + History_1 / 2.0f;
	o_Result = FlickerThere ? Interpolated : Current;
}		