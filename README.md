Please visit our docs [website](https://paulxshen.github.io/Luminescent.jl/) on running prebuilt binaries :) This repo now only hosts docs. 

# Home

Luminescent AI enables **generative design and simulation** of electromagnetic structures  in just a few lines of code! We help design next generation **photonic integrated circuits, optical metasurfaces, RF and microwave circuits, and antennas** in diverse industries including **consumer electronics, automotive, telecom, datacenters and quantum computing**. We created a **fully differentiable and GPU-accelerated finite difference time domain (FDTD)** simulator and geometry generator.

**Powerful**
- Generative inverse design and simulation in just few lines of Python code!
- Broadband and multimode S-parameters 
- Embedded mode solver for modal sources and monitors
- .gds and `gdsfactory` integration
- .stl / .step  3D geometry import

**Fast**
- GPU acceleration on NVIDIA, AMD, and Apple Silicon
- Adaptive graded mesh reduces cell count
- Tensor subpixel smoothing boosts accuracy

**Smart**
- Fully differentiable (native automatic differentiation in Julia)
- Simultaneous inverse design of multiple 2D and 3D structures
- Length scale controlled geometry optimizer with fabrication constraints

**Comprehensive** (some features require additional dev)
- Modal sources, plane waves, Gaussian beams, custom sources
- Oblique sources and monitors
- PML, periodic, Bloch, PEC boundaries
- Near and far field radiation patterns
- Nonlinear, dispersive and anisotropic materials_library
