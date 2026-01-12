# Inverse Design
## Overview
Workflow same as simulation. Only difference is passing additional args to `make` : `designs`, `targets`, `optimizer`. `solve` then runs optimization loop instead of single simulation. Each iteration has forward and backward pass. Forward pass is a simulation that has to save all intermediate fields to compute gradients in backward pass. Each iteration is >5x slower than single simulation. >40 iterations typically needed to converge. Altogether >200x slower than single simulation.
## Design regions
 ```{eval-rst}
    .. autofunction:: luminescent.opt.Design
```
## Targets
 ```{eval-rst}
    .. autofunction:: luminescent.opt.Target
```
## Optimizer
```{eval-rst}
    .. autofunction:: luminescent.opt.Optimizer
```