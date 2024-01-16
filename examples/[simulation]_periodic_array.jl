"""
simulation of plane wave on periodic array of dielectric pillars in front of metal of plane
"""

######################## 
# tunable configs #
###################
name = "periodic_array"
data_dir = "examples/utils"

"simulation params"
T = 20.0f0 # simulation duration in [periods]
nres = 32
Courant = 0.5f0 # Courant number

"recording params"
frameat = 1 / 16 # captures frame every _ periods
framerate = 16 # playback speed
########################

using UnPack, LinearAlgebra
include("../src/FDTDEngine.jl")
using .FDTDEngine
include("$data_dir/plot_recipes.jl")

F = Float32
name = "periodic_array"

dx = 1.0f0 / nres # pixel resolution in [wavelengths]
L = F[8, 8] # domain dimensions in [wavelengths]

# setup FDTD
polarization = :TMz
boundaries = [Periodic(2), PEC(1)] # unspecified boundaries default to PML
monitors = [Monitor(L / 2)]
# monitors = [Monitor([6, [3, 5]])]
sources = [PlaneWave(t -> exp(-((t - 1) / 1)^2) * cos(F(2π) * t), -1; Jz=1)]
fdtd_configs = setup(boundaries, sources, monitors, L, dx, polarization; F, Courant, T)
@unpack esz0, hsz0, dt, geometry_padding, field_padding, source_effects, monitor_instances, fields = fdtd_configs

ϵ1 = 1 #
ϵ2 = 2.25f0 # 
b = F.([norm([x, y] .- esz0 ./ 2) < 1 / dx for x = 1:esz0[1], y = 1:esz0[2]]) # circle
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

μ = ones(F, hsz0)
σ = zeros(F, esz0)
σm = zeros(F, hsz0)

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = [ϵ, μ, σ, σm]
u0 = collect(values(fields))

# run simulation
@showtime sol = accumulate((u, t) -> stepTMz(u, p, t, fdtd_configs), 0:dt:T, init=u0)
recordsim(sol, p, fdtd_configs, "$(name)_nres_$nres.mp4", title="$name"; frameat, framerate)
# plotmonitors(sol0, monitor_configs,)
