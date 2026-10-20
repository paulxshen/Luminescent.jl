# Usage
## Lumi Hub (cloud)
Run FDTD notebooks for free on [Lumi Hub Playground](https://hub.luminescentai.com/)! Similar to Google Colab but preloaded with Lumi FDTD, Copilot AI and persistent storage. Has runtime and size restrictions and excludes inverse design. Upgrade to Lumi Hub Enterprise for unrestricted usage: info@luminescentai.com
## Claude AI
Experimental feature letting Claude use Lumi FDTD! Say for example "design SOI220 ring resonator critically coupled Q>1000". Claude can optimize parameters and summarize results in a notebook. What can't Claude do?! First, [obtain token from Lumi Hub](https://hub.luminescentai.com/hub/token). Then, in Claude Code or Desktop (requires Claude Pro) connect our tool (replace YOUR_TOKEN)
`claude mcp add --transport http lumi https://api.luminescentai.com/mcp --header "Authorization: token YOUR_TOKEN"`

## Local machine
Recommended only for expert users. Backend still uses Lumi Hub server unless you bought binaries.
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
### Backend options
Lumi Hub server: `lumi.solve` sends jobs to Lumi Hub when `LUMI_API` is set. Use the token from [Lumi Hub](https://hub.luminescentai.com/hub/token). Lumi Hub's size and runtime limits apply.
```
export LUMI_API=https://api.luminescentai.com
export LUMI_TOKEN=YOUR_TOKEN
```
On Windows (PowerShell), then restart your terminal or IDE:
```
setx LUMI_API https://api.luminescentai.com
setx LUMI_TOKEN YOUR_TOKEN
```
Local binaries: executables can be compiled for enterprise customers above certain spend threshold. Buy once run forever no license server. Available for Linux and Windows.

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