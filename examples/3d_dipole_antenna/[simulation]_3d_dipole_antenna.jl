"""
simulation of dipole above ground plane
"""

using UnPack, LinearAlgebra, GLMakie
include("$(pwd())/src/main.jl")
include("$(pwd())/scripts/plot_recipes.jl")
# using FDTDEngine,FDTDToolkit


F = Float32
name = "3d_dipole"
T = 4.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]
Courant = 0.25f0 # Courant number

l = 2
sz0 = nres .* (l, l, l)

boundaries = [PEC(-3)] # unspecified boundaries default to PML
# n = [1, 0, 0]
monitors = []
sources = [
    Source(t -> t < 1 ? cos(F(2π) * t) : 0.0f0, [l / 2, l / 2, 0.25f0], [0, 0, 0.25f0]; Jz=1),
]
configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, ϵ, dt, geometry_padding, field_padding, source_effects, monitor_instances, fields, step, power = configs

p = apply(geometry_padding; ϵ, μ, σ, σm)
u0 = collect(values(fields))

# run simulation
@showtime sol = accumulate((u, t) -> step(u, p, t, configs), 0:dt:T, init=u0)
Ez = map(sol) do u
    u[3]
end
dir = @__DIR__
recordsim(Ez, p[1], configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=true)
