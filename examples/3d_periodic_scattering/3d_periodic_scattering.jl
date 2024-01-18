"""
simulation of plane wave scattering on periodic array 
"""

using UnPack, LinearAlgebra, GLMakie
include("$(pwd())/src/main.jl")
include("$(pwd())/scripts/plot_recipes.jl")
# using FDTDEngine,FDTDToolkit


F = Float32
name = "3d_scattering"
T = 4.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]
Courant = 0.25f0 # Courant number

l = 2
sz0 = nres .* (l, l, l)
ϵ1 = 1 #
ϵ2 = 2.25f0 # 
b = F.([norm(v .- sz0 ./ 2) < 0.5 / dx for v = Base.product(Base.oneto.(sz0)...)]) # sphere
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

boundaries = [Periodic(2), Periodic(3)]# unspecified boundaries default to PML
# n = [1, 0, 0]
monitors = []
sources = [
    PlaneWave(t -> t < 1 ? cos(F(2π) * t) : 0.0f0, -1; Jz=1)
]
configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, dt, geometry_padding, field_padding, source_effects, monitor_instances, fields, step, power = configs

p = apply(geometry_padding; ϵ, μ, σ, σm)
u0 = collect(values(fields))

# run simulation
@showtime sol = accumulate((u, t) -> step(u, p, t, configs), 0:dt:T, init=u0)
Ez = map(sol) do u
    u[3]
end
dir = @__DIR__
recordsim(Ez, p[1], configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=true)
