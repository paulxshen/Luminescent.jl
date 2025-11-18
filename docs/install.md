# Installation
Please fill out [new user form](https://forms.gle/fP9wAkdJinT8t66w8) to obtain free tier `__INSTALL_LINK__` which is pasted into scripts below. Free tier doesn't include inverse design.
## Google Colab Linux with CUDA
Easiest way is to run Luminescent on free Google Colab with GPU runtime. For more juice simply upgrade to Colab Pro to use A100 GPU for $1/hr. Run the following code cells.
```
import os
os.environ['LD_LIBRARY_PATH']='/temp-nvidia'
os.environ["JULIA_CUDA_USE_COMPAT"] = "0"
os.environ["PATH"] += ":/usr/local/Luminescent/bin"
```
```
%%shell
# frontend
apt-get install -y libglu1-mesa
pip install -U luminescent mediapy

# backend binaries
gdown --fuzzy __INSTALL_LINK__
tar -xf luminescent-latest.tar.gz  -C /usr/local/

# use prebundled CUDA but prevent LD_LIBRARY_PATH from interfering other packages
mkdir -p /temp-nvidia
mv /usr/lib64-nvidia/libcu* /temp-nvidia
mv /usr/lib64-nvidia/libnv* /temp-nvidia

# first time run downloads artifacts
luminescent
```
## Generic x86_64 Linux with CUDA
You can also run on your own machine or cloud provider. If install fails (probably due to CUDA config), raise an issue on [GitHub](https://github.com/paulxshen/Luminescent.jl).
```
# frontend
pip install -U luminescent 

# backend binaries
apt-get install -y libglu1-mesa
gdown --fuzzy __INSTALL_LINK__
tar -xf luminescent-latest.tar.gz  -C /usr/local/

export JULIA_CUDA_USE_COMPAT=0
export PATH=$PATH:/usr/local/Luminescent/bin

# first time run downloads artifacts
luminescent
```
## Windows, Apple, AMD ROCm
Please request