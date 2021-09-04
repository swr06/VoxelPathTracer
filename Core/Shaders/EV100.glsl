// implemented functions from the filament brdf docs.


#version 330 core

#define LENS_LENGTH 15.0
#define APERTURE 4.0
#define ISO 120.0
#define SHUTTER_SPEED 50.0
const float K = 12.5;
const float S = 100.0;

float computeEV100() {
     return log2((APERTURE * APERTURE) / (SHUTTER_SPEED) * 100.0f / (ISO));
}

float computeEV100fromLuma(float avgLuminance) {
     return log2(avgLuminance * (S / K));
}

float EV100ToExposure(float EV100) {
     return 1.0 / (exp2(EV100) * 1.2f);
}