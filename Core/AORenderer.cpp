#include "AORenderer.h"
#include <array>

float lerp(float a, float b, float f)
{
    return a + f * (b - a);
}

GLuint VoxelRT::SSAORenderer::GenerateSSAOKernelTexture()
{
    std::uniform_real_distribution<float> RandomGenerator(0.0, 1.0); 
    std::default_random_engine generator;
    std::array<glm::vec3, 64> SSAOKernel;

    SSAOKernel[0] = glm::vec3(0.04977, -0.04471, 0.04996);
    SSAOKernel[1] = glm::vec3(0.01457, 0.01653, 0.00224);
    SSAOKernel[2] = glm::vec3(-0.04065, -0.01937, 0.03193);
    SSAOKernel[3] = glm::vec3(0.01378, -0.09158, 0.04092);
    SSAOKernel[4] = glm::vec3(0.05599, 0.05979, 0.05766);
    SSAOKernel[5] = glm::vec3(0.09227, 0.04428, 0.01545);
    SSAOKernel[6] = glm::vec3(-0.00204, -0.0544, 0.06674);
    SSAOKernel[7] = glm::vec3(-0.00033, -0.00019, 0.00037);
    SSAOKernel[8] = glm::vec3(0.05004, -0.04665, 0.02538);
    SSAOKernel[9] = glm::vec3(0.03813, 0.0314, 0.03287);
    SSAOKernel[10] = glm::vec3(-0.03188, 0.02046, 0.02251);
    SSAOKernel[11] = glm::vec3(0.0557, -0.03697, 0.05449);
    SSAOKernel[12] = glm::vec3(0.05737, -0.02254, 0.07554);
    SSAOKernel[13] = glm::vec3(-0.01609, -0.00377, 0.05547);
    SSAOKernel[14] = glm::vec3(-0.02503, -0.02483, 0.02495);
    SSAOKernel[15] = glm::vec3(-0.03369, 0.02139, 0.0254);
    SSAOKernel[16] = glm::vec3(-0.01753, 0.01439, 0.00535);
    SSAOKernel[17] = glm::vec3(0.07336, 0.11205, 0.01101);
    SSAOKernel[18] = glm::vec3(-0.04406, -0.09028, 0.08368);
    SSAOKernel[19] = glm::vec3(-0.08328, -0.00168, 0.08499);
    SSAOKernel[20] = glm::vec3(-0.01041, -0.03287, 0.01927);
    SSAOKernel[21] = glm::vec3(0.00321, -0.00488, 0.00416);
    SSAOKernel[22] = glm::vec3(-0.00738, -0.06583, 0.0674);
    SSAOKernel[23] = glm::vec3(0.09414, -0.008, 0.14335);
    SSAOKernel[24] = glm::vec3(0.07683, 0.12697, 0.107);
    SSAOKernel[25] = glm::vec3(0.00039, 0.00045, 0.0003);
    SSAOKernel[26] = glm::vec3(-0.10479, 0.06544, 0.10174);
    SSAOKernel[27] = glm::vec3(-0.00445, -0.11964, 0.1619);
    SSAOKernel[28] = glm::vec3(-0.07455, 0.03445, 0.22414);
    SSAOKernel[29] = glm::vec3(-0.00276, 0.00308, 0.00292);
    SSAOKernel[30] = glm::vec3(-0.10851, 0.14234, 0.16644);
    SSAOKernel[31] = glm::vec3(0.04688, 0.10364, 0.05958);
    SSAOKernel[32] = glm::vec3(0.13457, -0.02251, 0.13051);
    SSAOKernel[33] = glm::vec3(-0.16449, -0.15564, 0.12454);
    SSAOKernel[34] = glm::vec3(-0.18767, -0.20883, 0.05777);
    SSAOKernel[35] = glm::vec3(-0.04372, 0.08693, 0.0748);
    SSAOKernel[36] = glm::vec3(-0.00256, -0.002, 0.00407);
    SSAOKernel[37] = glm::vec3(-0.0967, -0.18226, 0.29949);
    SSAOKernel[38] = glm::vec3(-0.22577, 0.31606, 0.08916);
    SSAOKernel[39] = glm::vec3(-0.02751, 0.28719, 0.31718);
    SSAOKernel[40] = glm::vec3(0.20722, -0.27084, 0.11013);
    SSAOKernel[41] = glm::vec3(0.0549, 0.10434, 0.32311);
    SSAOKernel[42] = glm::vec3(-0.13086, 0.11929, 0.28022);
    SSAOKernel[43] = glm::vec3(0.15404, -0.06537, 0.22984);
    SSAOKernel[44] = glm::vec3(0.05294, -0.22787, 0.14848);
    SSAOKernel[45] = glm::vec3(-0.18731, -0.04022, 0.01593);
    SSAOKernel[46] = glm::vec3(0.14184, 0.04716, 0.13485);
    SSAOKernel[47] = glm::vec3(-0.04427, 0.05562, 0.05586);
    SSAOKernel[48] = glm::vec3(-0.02358, -0.08097, 0.21913);
    SSAOKernel[49] = glm::vec3(-0.14215, 0.19807, 0.00519);
    SSAOKernel[50] = glm::vec3(0.15865, 0.23046, 0.04372);
    SSAOKernel[51] = glm::vec3(0.03004, 0.38183, 0.16383);
    SSAOKernel[52] = glm::vec3(0.08301, -0.30966, 0.06741);
    SSAOKernel[53] = glm::vec3(0.22695, -0.23535, 0.19367);
    SSAOKernel[54] = glm::vec3(0.38129, 0.33204, 0.52949);
    SSAOKernel[55] = glm::vec3(-0.55627, 0.29472, 0.3011);
    SSAOKernel[56] = glm::vec3(0.42449, 0.00565, 0.11758);
    SSAOKernel[57] = glm::vec3(0.3665, 0.00359, 0.0857);
    SSAOKernel[58] = glm::vec3(0.32902, 0.0309, 0.1785);
    SSAOKernel[59] = glm::vec3(-0.08294, 0.51285, 0.05656);
    SSAOKernel[60] = glm::vec3(0.86736, -0.00273, 0.10014);
    SSAOKernel[61] = glm::vec3(0.45574, -0.77201, 0.00384);
    SSAOKernel[62] = glm::vec3(0.41729, -0.15485, 0.46251);
    SSAOKernel[63] = glm::vec3(-0.44272, -0.67928, 0.1865);
    
    GLuint texture;

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, 8, 8, 0, GL_RGB, GL_FLOAT, &SSAOKernel[0]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    return texture;
}

GLuint VoxelRT::SSAORenderer::GenerateSSAONoiseTexture()
{
    std::uniform_real_distribution<float> RandomGenerator(0.0, 1.0);
    std::default_random_engine generator;
    std::vector<glm::vec2> SSAONoise;

    for (int i = 0; i < 64; ++i)
    {
        glm::vec2 sample(
            RandomGenerator(generator) * 2.0 - 1.0,
            RandomGenerator(generator) * 2.0 - 1.0
        );

        SSAONoise.push_back(sample);
    }

    GLuint texture;

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, 8, 8, 0, GL_RG, GL_FLOAT, &SSAONoise[0]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    return texture;
}
