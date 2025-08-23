# Luminescent AI - FDTD Simulation and Inverse Design

2025/08/23  
v1.1 
Paul Shen  
<pxshen@alumni.stanford.edu>  

# Summary

Luminescent AI enables generative design and simulation of electromagnetic structures  in just a few lines of code! We help design next generation photonic integrated circuits, optical metasurfaces, RF and microwave circuits, and antennas in diverse industries including consumer electronics, automotive, telecom, datacenters and quantum computing. We created an **automatic differentiation (AD) compatible** and **GPU-accelerated** finite difference time domain (FDTD) simulator and geometry generator.

Experimental release 🥼. Expect critters  🐛🐞

# Features
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
- Nonlinear, dispersive and anisotropic materials

# Examples
[Test drive any of our examples for free on Google Colab GPUs!](https://www.luminescentai.com/product) 
<!-- [Test drive any of our examples for free on Google Colab GPUs!](https://colab.research.google.com/drive/1NQ222-Odjz4Yg_ZguFyhLTMUplpKafgX?usp=sharing)  -->
<!-- # Product tiers

```python
``` -->

<!-- # RF and microwave circuits
## Simulation examples
![alt text](assets/sim-8.gif) ![alt text](assets/image-8.png) -->
<!-- ## Inverse design examples
### Microstrip patch antenna
### Microstrip bandpass filter
### 3D printed RF lens -->
