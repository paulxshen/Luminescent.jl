cd "C:\Users\pxshe\OneDrive\Desktop\New folder"
xcopy /Y makedocs\build\* docs
@REM xcopy /Y examples\*\*.mp4 makedocs\src\assets
xcopy /Y examples\periodic_scattering\*.mp4 makedocs\src\assets
xcopy /Y examples\quarter_wavelength_antenna\*.mp4 makedocs\src\assets
xcopy /Y examples\quarter_wavelength_antenna\*.png makedocs\src\assets
@REM xcopy /Y examples\slab_waveguide\*.mp4 makedocs\src\assets
xcopy /Y examples\inverse_design_waveguide_bend\*.mp4 makedocs\src\assets
