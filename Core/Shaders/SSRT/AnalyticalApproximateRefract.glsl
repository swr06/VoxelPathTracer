vec3 ApproximateRefractCoord(vec3 ViewSpacePosition, vec3 ViewSpacePositionUnderwater, vec3 ViewSpaceNormal, float ClipDepth) 
{
	float AirIdx = 1.00029f;
	float WaterIdx = 1.33333f;
	vec3 RefractedVector = normalize(refract(normalize(ViewSpacePosition), ViewSpaceNormal, AirIdx / WaterIdx));
	float RefractionAmount = distance(ViewSpacePosition, ViewSpacePositionUnderwater);
    vec3 ApproximateHitPosition = ViewSpacePosition + RefractedVector * RefractionAmount;
    vec3 HitCoordinate = ConvertToScreenSpace(ApproximateHitPosition);
    
	if(hitCoordinate.z < ClipDepth || HitCoordinate != clamp(HitCoordinate, 0.00001f, 0.99999f)) 
	{
            HitCoordinate = v_TexCoords;
    }
}