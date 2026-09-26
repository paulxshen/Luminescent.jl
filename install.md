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
sudo apt install -y harminv  libxi6 libxxf86vm1 libxrender1 libgl1 libglx-mesa0 libsm6 libice6 libxext6 libglib2.0-0 libglu1-mesa libxcursor1 libxft2 libxinerama1
```
### Linux ARM
request
### Windows
Prerequisite: [Microsoft Visual C++ Redistributable](https://aka.ms/vc14/vc_redist.x64.exe)
Powershell (run as Administrator):
```
$ProgressPreference = 'SilentlyContinue'
Invoke-WebRequest -Uri "https://storage.googleapis.com/lumi-fdtd/lumi.zip" -OutFile "$env:TEMP\lumi.zip"
Expand-Archive -Path "$env:TEMP\lumi.zip" -DestinationPath "C:\Program Files" -Force

[Environment]::SetEnvironmentVariable(
    "Path",
    [Environment]::GetEnvironmentVariable("Path", "User") + ";C:\Program Files\Luminescent\bin",
    "User"
)
```
### MacOS
request
 
<!-- - Backend: install Julia. Then `julia -e 'using Pkg; Pkg.add(url="https://github.com/paulxshen/Luminescent.jl")'` 

If GPU, also install CUDA.jl: `julia -e 'using Pkg; Pkg.add("CUDA")'` and use `lm.solve(..., backend="cuda")` in Python -->