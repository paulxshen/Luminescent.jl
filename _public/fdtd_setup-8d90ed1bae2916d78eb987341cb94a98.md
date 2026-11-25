# Install
| Tier | Simulation | Inverse design | Compute |
|---|---|---|---|
| Free | yes | no | fractional GPU |
| Pro | yes | no | full GPU |
| Max | yes | yes | multi-GPU |

Upgrade to Pro or Max: info@luminescentai.com


## Option 1: Local Python + Remote Lumi Colab 
Run Python locally. `lumi.solve` uploads problem files to Lumi Colab server and downloads solution files. 
- `pip install -U luminescent`
- Get your [Lumi Colab token](https://colab.luminescentai.com/hub/token), then save it once: `python -c "import luminescent as lumi; lumi.set_token('<token>')"`
## Option 2: fully local 
Local Python and FDTD executable. Available to paid customers above certain spend threshold. 
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