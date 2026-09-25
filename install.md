# Installation
Install both frontend Python API and backend FDTD executable.
The free version has most simulation features but has size limit and excludes inverse design.
## Frontend
 ```
 pip install luminescent
 ``` 
## Backend
### Linux x86_64
```
wget https://storage.googleapis.com/lumi-fdtd/lumi.tar.gz -O /tmp/lumi.tar.gz 
sudo tar -x -f /tmp/lumi.tar.gz -C /usr/local --strip-components 1
sudo apt-get update
sudo apt install harminv

# Optional: in case running on cloud VM or headless server, you might need:
sudo apt install -y libxi6 libxxf86vm1 libxrender1 libgl1 libglx-mesa0 libsm6 libice6 libxext6 libglib2.0-0 libglu1-mesa libxcursor1 libxft2 libxinerama1
```
### Linux ARM
request
### Windows
Powershell (run as Administrator):
```
Invoke-WebRequest -Uri "https://storage.googleapis.com/lumi-fdtd/lumi.zip" -OutFile "$env:TEMP\lumi.zip"
Expand-Archive -Path "$env:TEMP\lumi.zip" -DestinationPath "C:\Program Files\lumi" -Force

[Environment]::SetEnvironmentVariable(
    "Path",
    [Environment]::GetEnvironmentVariable("Path", "User") + ";C:\Program Files\lumi\lumi",
    "User"
)
```
### MacOS
request
 
<!-- - Backend: install Julia. Then `julia -e 'using Pkg; Pkg.add(url="https://github.com/paulxshen/Luminescent.jl")'` 

If GPU, also install CUDA.jl: `julia -e 'using Pkg; Pkg.add("CUDA")'` and use `lm.solve(..., backend="cuda")` in Python -->