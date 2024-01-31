"""
simulation of quarter wavelength antenna above ground plane
"""

using UnPack, LinearAlgebra, GLMakie
# using FDTDEngine,FDTDToolkit
dir = pwd()
include("$(dir)/src/main.jl")
include("$dir/../FDTDToolkit.jl/src/main.jl")




F = Float32
name = "3d_quarter_wavelength_antenna"
T = 8.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]
Courant = 0.8 / √3 # Courant number

l = 2
sz0 = nres .* (l, l, l)

boundaries = [] # unspecified boundaries default to PML
# boundaries = [PEC(-3)] # unspecified boundaries default to PML
# n = [1, 0, 0]
monitors = []
sources = [
    Source(t -> cos(F(2π) * t), [l / 2, l / 2, 0.125f0], [0, 0, 0.25f0]; Jz=1),
]
configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, ϵ, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, fields, power = configs
# volume(source_effects[1]._g[:Jz])

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)
u0 = collect(values(fields))

# run simulation
t = 0:dt:T
sol = similar([u0], length(t))
@showtime sol = accumulate!((u, t) -> step!(u, p, t, configs), sol, t, init=u0)
Ez = map(sol) do u
    u[3]
end
dir = @__DIR__
@showtime recordsim(Ez, p[1][3], configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=true)