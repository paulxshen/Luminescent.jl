# Features
**Fully differentiable and GPU-accelerated finite difference time domain (FDTD)** engine for full wave electromagnetic simulation and inverse design.

**Powerful**
- Generative inverse design and simulation in just few lines of Python code!
- Broadband and multimode S-parameters 
- Embedded mode solver for modal sources and monitors
- .gds and `gdsfactory` integration
- .stl / .step  3D geometry import

**Smart**
- Fully differentiable (native automatic differentiation in Julia)
- Simultaneous inverse design of multiple 2D and 3D structures
- Length scale controlled geometry optimizer with fabrication constraints

**Fast**
- GPU acceleration on NVIDIA, AMD, and Apple Silicon
- Adaptive graded mesh reduces cell count
- Tensor subpixel smoothing boosts accuracy

**Comprehensive**
- Modal sources, plane waves, Gaussian beams, custom sources
- Oblique sources and monitors
- PML, periodic, PEC, PMC boundaries
- Near and far field radiation patterns
- Nonlinear, dispersive and anisotropic materials

![](simulation.gif)

Fig. Inverse designed broadband bidirectional perfectly vertical grating coupler (PVGC) 