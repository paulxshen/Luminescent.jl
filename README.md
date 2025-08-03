# Luminescent AI - FDTD Simulation and Inverse Design

2025/08/01  
Paul Shen  
<pxshen@alumni.stanford.edu>  

# Summary

Luminescent AI enables generative design and simulation of electromagnetic structures  in just a few lines of code! We help design next generation photonic integrated circuits, optical metasurfaces, RF and microwave circuits, and antennas in diverse industries including consumer electronics, automotive, telecom, datacenters and quantum computing. We created an automatic differentiation (AD) and GPU compatible finite difference time domain (FDTD) simulator and geometry generator.


[GitHub](https://github.com/paulxshen/Luminescent.jl): Star us :) We respond to issues within a day  
[LinkedIn](https://www.linkedin.com/company/luminescent-ai/about): Follow us for new features and bug fixes  
[Company](luminescentai.com): Consulting, collaboration, publication, investment

Email: pxshen@alumni.stanford.edu  
WhatsApp: 650-776-7724  
WeChat: pxshen1230  

Experimental release 🥼. Expect critters  🐛🐞

# Features
 
**Easy**
- Generative design and simulation in just few lines of Python code!
- Broadband and multimode S-parameters 
- Embedded mode solver for modal sources and monitors
- Near and far field radiation patterns
- Application specific API via gdsfactory integration
    - Photonic integrated circuits 
    - Optical metasurfaces
    - RF and microwave circuits: microstrip, CPWG
    - Planar antennas
- General API via .stl / .step  3D geometry import
    - 3D antennas
    - Biomedical RF

**Fast**
- Adaptive graded mesh for reduced cell count
- Tensor subpixel smoothing for accuracy even on coarse mesh
- Optional GPU acceleration

**Smart**
- Generative AI design of 2D and 3D structures
- Fully differentiable for inverse design
- Automatic differentiation (AD) gradients for adjoint optimization
- Length scale controlled geometry optimizer with fabrication constraints  

**Comprehensive**
- Modal sources, plane waves, Gaussian beams, custom sources
- Oblique sources and monitors
- PML, periodic, Bloch, PEC boundaries
- Nonlinear, dispersive and anisotropic materials

# Examples
## Simulation
### [Tutorial: photonic ring resonator](https://colab.research.google.com/drive/1NQ222-Odjz4Yg_ZguFyhLTMUplpKafgX?usp=sharing)  
### [Microwave frequency selective surface (FSS)](https://colab.research.google.com/drive/1R477Sk5y-qPRHokCeYm98U1Jo4C8qwR8?usp=sharing) | [video](https://www.youtube.com/watch?v=Uq7OnLmFSEk)
### [Microstrip quarter wavelength stub filter](https://colab.research.google.com/drive/1NqAf7z9sDB8BdkNkSDoP24ckSv0XYquH?usp=sharing)

## Inverse design 
### [Tutorial: photonic wavelength demultiplexer](https://colab.research.google.com/drive/1vJ-y2NabANt9xD-hWUq0rTBaeLzaBzmG?usp=sharing) | [video](https://youtu.be/DX1_vPu0ngI)

### [Photonic symmetric crossing](https://colab.research.google.com/drive/1ah2517VuLNj0kUSyux2LGEtQppVClvNH?usp=sharing)

# Product
All installers are fully local and never expire.   

**Free**
- CPU + GPU binaries
- 3D simulation up to 0.5M cells
- 2.5D inverse design  

**Pro**
- Everything in free
- Unrestricted 3D simulation
- 3D inverse design
- 1 year of updates and support
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
