"""
Solves first photonics inverse design problem in Google's Ceviche challenges. 
Design is a tight right angle waveguide bend.
"""

######################## 
# tunable configs #
###################
name = "waveguide_bend"
dir = "examples/utils"

"training params"
T = 4.0f0 # simulation duration in [periods]
nres = 32
dx = 1.0f0 / nres
Courant = 0.25f0 # Courant number

########################

using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Flux, Optim, Jello

include("../src/main.jl")
# include("../src/FDTDEngine.jl")
# using .FDTDEngine

include("../../FDTDToolkit.jl/src/plot_recipes.jl")

F = Float32
l = 1
sz0 = l .* (nres, nres, nres)
ϵ = ones(F, sz0)

boundaries = [] # unspecified boundaries default to PML
boundaries = [Periodic(2), Periodic(3)]
# n = [1, 0, 0]
monitors = []
sources = [
    # UniformSource(t -> cos(F(2π) * t), l/2*[1,1,1], [0, 0, 0.25f0]; Jz=1),
    PlaneWave(t -> exp(-((t - 1) / 1)^2) * cos(F(2π) * t), -1; Jz=1)
]
fdtd_configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, dt, geometry_padding, field_padding, source_effects, monitor_instances, fields = fdtd_configs




p = apply(geometry_padding; ϵ, μ, σ, σm)
u0 = collect(values(fields))

# run simulation
@showtime sol = accumulate((u, t) -> step3(u, p, t, fdtd_configs), 0:dt:T, init=u0)
nres = 1 ÷ dx
# E = map(sol) do u
#     u[1] .^ 2 + u[2] .^ 2 + u[3] .^ 2
# end
# recordsim(E, nothing, fdtd_configs, "$(name)_nres_$nres.mp4",umax=.5, title="$name"; playback=.25, bipolar=false)
E = map(sol) do u
    u[3]
end
recordsim(E, nothing, fdtd_configs, "$(name)_nres_$nres.mp4", title="$name"; playback=0.25, bipolar=true)
# plotmonitors(sol0, monitor_configs,)
display(plotstep(E[16], nothing, fdtd_configs))
display(volume(E[32]))
