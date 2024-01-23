"""
simulation of dipole above ground plane
"""

using UnPack, LinearAlgebra, GLMakie
include("$(pwd())/src/main.jl")
include("$(pwd())/scripts/plot_recipes.jl")
# using FDTDEngine,FDTDToolkit


F = Float32
name = "3d_dipole"
T = 2.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]
Courant = 0.25f0 # Courant number

l = 2
sz0 = nres .* (l, l, l)

boundaries = [] # unspecified boundaries default to PML
# boundaries = [PEC(-3)] # unspecified boundaries default to PML
# n = [1, 0, 0]
monitors = []
sources = [
    Source(t -> t < 1 ? cos(F(2π) * t) : 0.0f0, [l / 2, l / 2, l / 2], [0, 0, 0.25f0]; Jz=1),
]
configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, ϵ, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, fields, step, power = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)
u0 = collect(values(fields))

# run simulation
@showtime sol = accumulate((u, t) -> step(u, p, t, configs), 0:dt:T, init=u0)
Ez = map(sol) do u
    u[3]
end
dir = @__DIR__
recordsim(Ez, p[1][3], configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=true)
volume(Ez[16])
volume(p[3][3])
volume(σ)
volume(source_effects[1]._g[:Jz])
Ez = ones(16, 16, 16)
Ez, = apply(field_padding; Ez)
volume(Ez)