#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

// Packed Texture, xyz contains color and the w component contains depth
uniform sampler2D u_Texture;

uniform vec2 u_Dimensions;

const vec2 BokehOffsets[50] = vec2[50] ( vec2(1.5017174006515064e-17, 0.24524906310898137), vec2(-0.2332457195850781, 0.07578612835520925), vec2(-0.1441537824340053, -0.1984106599096999), vec2(0.14415378243400523, -0.19841065990969994), vec2(0.23324571958507806, 0.07578612835520915), vec2(1.5017174006515063e-16, 0.49049812621796274), vec2(-0.25884016390759723, 0.35626292189421216), vec2(-0.4664914391701561, 0.15157225671041863), vec2(-0.41881218285608623, -0.13608032723223978), vec2(-0.2883075648680107, -0.39682131981939966), vec2(-1.887521368465203e-16, -0.4403651893239443), vec2(0.28830756486801035, -0.3968213198193999), vec2(0.4188121828560862, -0.13608032723224017), vec2(0.46649143917015623, 0.15157225671041819), vec2(0.25884016390759684, 0.35626292189421244), vec2(1.7124132500063345e-15, 0.7357471893269439), vec2(-0.27506591574761147, 0.6178081620167741), vec2(-0.4908826676496936, 0.4419927392274961), vec2(-0.6997371587552336, 0.22735838506562794), vec2(-0.672570520845704, -0.07069001025572542), vec2(-0.5720511613453199, -0.33027389199295754), vec2(-0.4324613473020151, -0.5952319797291), vec2(-0.1406055259662536, -0.661496991019891), vec2(0.13733560663443514, -0.6461132300758565), vec2(0.43246134730201663, -0.5952319797290995), vec2(0.5856715267925036, -0.33813761365035117), vec2(0.6569292341109864, -0.06904604477490263), vec2(0.6997371587552346, 0.22735838506562583), vec2(0.5025704357670656, 0.45251645290919346), vec2(0.26866898824959135, 0.6034404276162195), vec2(-9.617063541017877e-16, 0.9809962524359257), vec2(-0.29009284723614176, 0.8928139801909716), vec2(-0.5176803278151952, 0.7125258437884238), vec2(-0.7154325033838056, 0.5197921396256749), vec2(-0.9329828783403125, 0.303144513420836), vec2(-0.9387601734426196, 8.047547491014627e-16), vec2(-0.8376243657121732, -0.2721606544644777), vec2(-0.7154325033838066, -0.5197921396256738), vec2(-0.5766151297360206, -0.7936426396388002), vec2(-0.29009284723614365, -0.8928139801909724), vec2(7.555537024365194e-16, -0.8807303786478886), vec2(0.27327089963618906, -0.8410413489993727), vec2(0.5766151297360188, -0.793642639638801), vec2(0.7594729339574509, -0.5517893853890959), vec2(0.8376243657121717, -0.27216065446448223), vec2(0.8843232074952327, -8.663868683766059e-16), vec2(0.9329828783403126, 0.3031445134208341), vec2(0.7594729339574525, 0.5517893853890949), vec2(0.5176803278151966, 0.7125258437884228), vec2(0.27327089963619067, 0.8410413489993722) );
float BLUR_SIZE = 20.0f;
const float BLUR_SCALE = 10.0f;
const float VIRTUAL_FAR_PLANE = 300.0f;

float GetBlurRadius(float SampleDepth, float FocusPoint, float FocusScale)
{
	float CircleOfConfusion = clamp((1.0f / FocusPoint - 1.0f / SampleDepth) * FocusScale, -1.0f, 1.0f);
    return abs(CircleOfConfusion) * BLUR_SIZE;
}

float GetCircleOfConfusion(float Depth, float CenterDepth, float Scale) 
{
	float CircleOfConfusion = abs(Depth - CenterDepth) / 0.6;
	return CircleOfConfusion / (1.0f / Scale + CircleOfConfusion);
}

void main() {
	vec2 TexelSize = 1.0f / textureSize(u_Texture, 0).xy;

	vec4 CenterSample = texture(u_Texture, v_TexCoords);
	float CenterDepth = 10.0f / CenterSample.w;

	vec3 TotalColor = CenterSample.xyz;
	float TotalWeight = 1.0f;

	float FocusPoint = 10.0f / texture(u_Texture, vec2(0.5f)).w;

	float CenterRadius = GetBlurRadius(CenterDepth, FocusPoint, BLUR_SCALE);

	float RadiusScale = 1.0f;

	for (int Sample = 0 ; Sample < 50 ; Sample++) {
		
		vec2 Offset = BokehOffsets[Sample];
		vec2 SampleCoord = v_TexCoords + Offset * TexelSize;
		vec4 SampleData = texture(u_Texture, SampleCoord);
		float SampleDepth = 10.0f / SampleData.w;

		float BlurRadiusAt = GetBlurRadius(SampleDepth, FocusPoint, BLUR_SCALE);

		if (SampleDepth > CenterDepth) {
			BlurRadiusAt = clamp(BlurRadiusAt, 0.0f, CenterRadius * 2.0f);
		}	

		float m = smoothstep(RadiusScale - 0.5f, RadiusScale + 0.5f, BlurRadiusAt);

		TotalColor += mix(TotalColor / TotalWeight, SampleData.xyz, m);
		TotalWeight += 1.0f;
	}

	o_Color = TotalColor / TotalWeight;
}