magick.exe convert %CD%/BASEPBR.png -channel RGBA -separate SEPARATED.png
ren "%CD%\SEPARATED-0.png" "Roughness.png"
ren "%CD%\SEPARATED-1.png" "Metalness.png"
magick.exe convert %CD%/Roughness.png -channel R -negate Roughness.png
del %CD%\SEPARATED-2.png
del %CD%\SEPARATED-3.png
magick.exe convert %CD%/Normal.png -channel B -separate Displacement.png
magick.exe convert %CD%/Roughness.png %CD%/Metalness.png %CD%/Displacement.png %CD%/AO.png -channel RGBA -combine PBR.png