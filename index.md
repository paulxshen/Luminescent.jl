# Home
## Overview
Lumi FDTD (finite difference time domain) powers full wave electromagnetic simulation and inverse design in photonics and RF. Maintained by Luminescent AI founded by Stanford alum designing next generation photonic chips. Design consulting also available. Below example is AI designed perfectly vertical grating coupler:

![](simulation.gif)  


Raise Issues at [GitHub](https://github.com/paulxshen/Luminescent.jl)  
Follow on [LinkedIn](https://www.linkedin.com/company/luminescent-ai)  
Email: pxshen@alumni.stanford.edu  
WhatsApp and WeChat: +1 (650) 776-7724

# Features

**AI design**  
- Fully differentiable for inverse design
- Topology optimization of multiple 2D or 3D structures
- Length scale constrained for fabrication compliance

**Comprehensive**
- Anisotropic, nonlinear, dispersive, metallic materials
- PML, periodic, PEC, PMC boundaries
- Adaptive mesh for thin or sharp geoemtries 
- Modal or field monitors in time or frequency domain
- Embedded FDFD and FEM mode solver 
<!-- - Near and far field radiation patterns -->

**Easy**
- Inverse design and simulation in Python frontend 
- .gds KLayout and `gdsfactory` integration
- .stl / .step 3D geometry import

**Performant**
- Adaptive mesh for cell count efficiency
- Tensor subpixel smoothing for accuracy
- Optional GPU acceleration on NVIDIA, AMD or Apple

**Coming soon**
- Bloch boundary
- RF features
-- Radiation patterns
-- farfield projections