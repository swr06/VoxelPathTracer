// Implementation of screen space contact shadows
// Might make this stochastic to get nice soft shadows 

#version 430 core

in vec2 v_TexCoords; // Screen space uv

layout (location = 0) out float o_Shadow; 

// Depth
uniform sampler2D u_DepthTexture;

// Matrices 
uniform mat4 u_ProjectionMatrix;
uniform mat4 u_InverseProjectionMatrix;
uniform mat4 u_ViewMatrix;
uniform mat4 u_InverseViewMatrix;

// Near and far plane 
uniform float u_zNear;
uniform float u_zFar;

// Trace direction, ie the sun direction
uniform vec3 u_SunDir;

vec3 ViewPosFromDepth(float depth)
{
    float z = depth * 2.0f - 1.0f; 

    vec4 ClipSpacePosition = vec4(v_TexCoords * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjectionMatrix * ClipSpacePosition;

    ViewSpacePosition /= ViewSpacePosition.w;

    return ViewSpacePosition.xyz;
}

//Z buffer is nonlinear by default, so we fix this here
float linearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

vec3 ViewSpaceToClipSpace(in vec3 view_space)
{
	vec4 clipSpace = u_ProjectionMatrix * vec4(view_space, 1);
	vec3 NDCSpace = clipSpace.xyz / clipSpace.w;
	vec3 screenSpace = 0.5 * NDCSpace + 0.5;
	return screenSpace;
}

vec3 ViewSpaceToScreenSpace(vec3 ViewSpace) 
{
    vec4 ClipSpace = u_ProjectionMatrix * vec4(ViewSpace, 1.0);
    ClipSpace.xyz = ClipSpace.xyz / ClipSpace.w;
    return ClipSpace.xyz * 0.5f + 0.5f;
}

vec3 convertScreenSpaceToWorldSpace(in vec2 txc)
{
    float d = texture(u_DepthTexture, txc).r;
    return ViewPosFromDepth(d);
}

bool WithinBounds(in vec3 p)
{
    if (p.x < 0.0f || p.x > 1.0f ||
		p.y < 0.0f || p.y > 1.0f ||
		p.z < 0.0f || p.z > 1.0f)
	{
		return false;
	}

    return true;
}


float RaytraceShadow()
{
    const int MAX_REFINEMENTS = 3;
    float InitialStepAmount = 1.0 - clamp(0.1f / 100.0, 0.0, 0.99);
    vec3 ViewSpacePosition = ViewPosFromDepth(texture(u_DepthTexture, v_TexCoords).r);
    vec3 ViewSpaceViewDirection = normalize(ViewSpacePosition);
    vec3 NormalizedSunDir = normalize(u_SunDir);
    vec3 ViewSpaceVector = InitialStepAmount * NormalizedSunDir;
    vec3 ViewSpaceVectorFAR = u_zFar * NormalizedSunDir; // Scaled by far plane
	vec3 PreviousPosition = ViewSpacePosition;
    vec3 ViewSpaceVectorPosition = PreviousPosition + ViewSpaceVector;
    vec3 CurrentPosition = ViewSpaceToScreenSpace(ViewSpaceVectorPosition);
	int NumRefinements = 0;
	vec2 FinalUV = vec2(-1.0f);

	float finalSampleDepth = 0.0;

    // march!
    for (int Step = 0; Step < 40; Step++)
    {
        if(-ViewSpaceVectorPosition.z > u_zFar * 1.4f || -ViewSpaceVectorPosition.z < 0.0f)
        {
		    break;
		}

        vec2 SamplePos = CurrentPosition.xy;
        float SampleDepth = convertScreenSpaceToWorldSpace(SamplePos).z;
        float CurrentDepth = ViewSpaceVectorPosition.z;
        float diff = SampleDepth - CurrentDepth;
        float error = length(ViewSpaceVector / pow(2.0f, NumRefinements));

        if(diff >= 0 && diff <= error * 2.0f && NumRefinements <= MAX_REFINEMENTS)
        {
        	ViewSpaceVectorPosition -= ViewSpaceVector / pow(2.0f, NumRefinements);
        	NumRefinements++;
		}

		else if (diff >= 0 && diff <= error * 4.0f && NumRefinements > MAX_REFINEMENTS)
		{
			FinalUV = SamplePos;
			finalSampleDepth = SampleDepth;
			break;
		}

        ViewSpaceVectorPosition += ViewSpaceVector / pow(2.0f, NumRefinements);

        if (Step > 1)
        {
            ViewSpaceVector *= 1.375f; // exponential step, works super well.
        }

		CurrentPosition = ViewSpaceToScreenSpace(ViewSpaceVectorPosition);

		if (!WithinBounds(CurrentPosition)) 
        {
            break;
        }
    }

    return FinalUV == clamp(FinalUV, 0.00001, 0.99999f) ? 1.0f : 0.0f;
}

void main()
{
    o_Shadow = RaytraceShadow();
}

//