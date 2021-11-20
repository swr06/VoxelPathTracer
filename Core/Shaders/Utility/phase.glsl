float rayleighPhase(float cosTheta) 
{
	float y = 0.035 / (2.0 - 0.035);
	float p1 = 3.0 / (4.0 * (1.0 + 2.0*y));
	float p2 = (1.0 + 3.0*y) + (1.0 - y) * square(cosTheta);
	float phase = p1 * p2;
	phase *= rcp(pi*4.0);
	return phase;
}

float cornetteShanksMiePhase(float cosTheta, float g) 
{
	float gg = g*g;
	float p1 = 3.0 * (1.0 - gg) * rcp((pi * (2.0 + gg)));
	float p2 = (1.0 + square(cosTheta)) * rcp(pow((1.0 + gg - 2.0 * g * cosTheta), 3.0/2.0));
	float phase = p1 * p2;
	phase *= rcp(pi*4.0);
	return max(phase, 0.0);
}
