##

Engineers run simulations to improve designs. Each time the design changes, the simulation is re-run. This can be done systematically in "parameter sweeps" where different combinations of parameter values are simulated to determine the best design. However, this scales exponentially wrt the number of parameters or DOFs. 

## General workflow 
We use gradient descent, the same as in machine learning. In lieu of optimizing neural network parameters, we're optimizing geometry (or source) parameters. In each training iteration, we generate geometry, run the simulation, calculate the objective metric, and do a backward pass to derive the gradient wrt the geometry parameters. We then do a gradient based parameter update in preparation for the next iteration.

The geometry is thus the first update. It typically has a static component which we can't change such as interfacing waveguides. Then there's a design component which we can change or optimize. The user is responsible for generating the design geometry wrt design parameters. If any pattern is allowed in the design region, our sister package `Jello.jl` can be used as a length scale controlled geometry generator. In any case, the result needs to be a 2d/3d array of each relevant materials property eg permitivity. 

With geometry ready, we can run the simulation. Duration is roughly the time it takes to reach steady state, such as how long it take for the signal to reach output port. The objective is usually a steady state metric which can be computed using values from the final period. 
We optimize geometry for some objective. 