#version 330 core

in vec2 v_TexCoords;

//Output
out vec3 o_Color;

//Uniforms
uniform sampler2D u_NormalTexture;
uniform sampler2D u_DepthTexture;
uniform sampler2D u_RefractionMaskTexture;

//Projection matrix
uniform mat4 u_ProjectionMatrix;
uniform mat4 u_InverseProjectionMatrix;
uniform mat4 u_ViewMatrix;
uniform float u_zNear;
uniform float u_zFar;

//Tweakable variables
const float InitialStepAmount = 0.025f;
const float StepRefinementAmount = 0.26f;
const int MaxRefinements = 16;

vec3 ViewPosFromDepth(float depth)
{
    float z = depth * 2.0f - 1.0f; // No need to linearize

    vec4 ClipSpacePosition = vec4(v_TexCoords * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjectionMatrix * ClipSpacePosition;

    // Perspective division
    ViewSpacePosition /= ViewSpacePosition.w;

    return ViewSpacePosition.xyz;
}

//Z buffer is nonlinear by default, so we fix this here
float linearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

vec3 ViewSpaceToScreenSpace(in vec3 view_space)
{
	vec4 clipSpace = u_ProjectionMatrix * vec4(view_space, 1);
	vec3 NDCSpace = clipSpace.xyz / clipSpace.w;
	vec3 screenSpace = 0.5 * NDCSpace + 0.5;
	return screenSpace;
}

vec2 ComputeRefraction()
{	
	//Values from textures
	vec2 ScreenSpacePosition2D = v_TexCoords;
	vec3 ViewSpacePosition = ViewPosFromDepth(texture(u_DepthTexture, ScreenSpacePosition2D).r);
	vec3 ViewSpaceNormal = vec3(u_ViewMatrix * vec4(texture(u_NormalTexture, ScreenSpacePosition2D).xyz, 0.0f));

	//Screen space vector
	vec3 ViewSpaceVector = (refract(normalize(ViewSpacePosition), normalize(ViewSpaceNormal), 1.0f / 5.73f));
	vec3 ScreenSpacePosition = ViewSpaceToScreenSpace(ViewSpacePosition);
	ScreenSpacePosition.z = linearizeDepth(ScreenSpacePosition.z);

	vec3 ViewSpaceVectorPosition = ViewSpacePosition + ViewSpaceVector;
	vec3 ScreenSpaceVectorPosition = ViewSpaceToScreenSpace(ViewSpaceVectorPosition);
	ScreenSpaceVectorPosition.z = linearizeDepth(ScreenSpaceVectorPosition.z);

	vec3 ScreenSpaceVector = InitialStepAmount * normalize(ScreenSpaceVectorPosition - ScreenSpacePosition);
	
	vec3 OldPosition = ScreenSpacePosition + ScreenSpaceVector;
	vec3 CurrentPosition = OldPosition + ScreenSpaceVector;
	
	//State
	vec4 color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	vec2 final_uv = vec2(-1.0f);
	int count = 0;
	int NumRefinements = 0;
	int depth = 0;

	for (int count = 0 ; count < 50 ; count++)
	{
		if(CurrentPosition.x < 0 || CurrentPosition.x > 1 ||
			CurrentPosition.y < 0 || CurrentPosition.y > 1 ||
			CurrentPosition.z < 0 || CurrentPosition.z > 1)
		{
			break;
		}

		//intersections
		vec2 SamplePos = CurrentPosition.xy;
		float CurrentDepth = (CurrentPosition.z);
		float SampleDepth = linearizeDepth(texture(u_DepthTexture, SamplePos).x);
		float diff = CurrentDepth - SampleDepth;

		if(diff >= 0 && diff < 1000.0f)
		{
			ScreenSpaceVector *= StepRefinementAmount;
			CurrentPosition = OldPosition;
			NumRefinements++;

			if(NumRefinements >= MaxRefinements)
			{
				final_uv = SamplePos;
				break;
			}
		}

		//Step ray
		OldPosition = CurrentPosition;
		CurrentPosition = OldPosition + ScreenSpaceVector;
	}

	return final_uv;
}

void main()
{
	o_Color = vec3(-1.0f);

	if (texture(u_RefractionMaskTexture, v_TexCoords).r == 1.0f)
	{
		vec2 RefractedUV = ComputeRefraction();
		o_Color = vec3(RefractedUV, 0.0f);
	}
}
