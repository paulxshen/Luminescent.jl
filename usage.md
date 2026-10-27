# FDTD Install
## Option 1: Lumi Colab 
Run FDTD notebooks for free on [Lumi Colab](https://hub.luminescentai.com/)! Similar to Google Colab but preloaded with Lumi FDTD, Copilot AI and persistent storage. Has runtime and size restrictions and excludes inverse design. Upgrade to Lumi Colab Enterprise tier for unrestricted usage: info@luminescentai.com


## Option 2: Local Python frontend + Lumi Colab backend
Run Python notebooks locally. `lumi.solve` uploads problem files to Lumi Colab server and downloads solution files. Same tiers apply
### Prereqs
- Python (Anaconda distro recommended for Windows users)
- VS Code or another IDE
- On a Linux VM you may need
```
sudo apt-get update
sudo apt install -y harminv  libxi6 libxxf86vm1 libxrender1 libgl1 libglx-mesa0 libsm6 libice6 libxext6 libglib2.0-0 libglu1-mesa libxcursor1 libxft2 libxinerama1
```
### Frontend
`pip install luminescent`
### Backend
Set environmental variables with YOUR_TOKEN from [Lumi Colab](https://hub.luminescentai.com/hub/token). A copy of each run is kept in your Lumi Colab home at `~/local/runs/<run folder name>`.
On Linux:  
```
export LUMI_API=https://api.luminescentai.com
export LUMI_TOKEN=YOUR_TOKEN
```
On Windows (PowerShell), then restart your terminal or IDE:
```
setx LUMI_API https://api.luminescentai.com
setx LUMI_TOKEN YOUR_TOKEN
```
## Option 3: fully local 
Executables can be compiled for enterprise customers above major spend threshold. 
<!-- Only relevant for customers who  -->
<!-- Buy once run forever no license server.  -->

<!-- sudo python -e 'import luminescent' # installs backend for first time -->
<!-- Install Windwos Subsystem for Linux (WSL) Ubuntu. In VS Code cd to Ubuntu drive and follow previous Linux backend instructions.  -->
<!-- Prerequisite: [Microsoft Visual C++ Redistributable](https://aka.ms/vc14/vc_redist.x64.exe)
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
``` -->
 
<!-- - Backend: install Julia. Then `julia -e 'using Pkg; Pkg.add(url="https://github.com/paulxshen/Luminescent.jl")'` 

If GPU, also install CUDA.jl: `julia -e 'using Pkg; Pkg.add("CUDA")'` and use `lm.solve(..., backend="cuda")` in Python -->