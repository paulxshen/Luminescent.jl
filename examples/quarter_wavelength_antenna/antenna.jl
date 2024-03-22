"""
simulation of quarter wavelength antenna above conductor ground plane
"""

using UnPack, LinearAlgebra, GLMakie
using Luminescent, LuminescentVisualization

# dir = pwd()
# include("$(dir)/src/main.jl")
# include("$dir/../LuminescentVisualization.jl/src/main.jl")

name = "quarter_wavelength_antenna"
F = Float32
dogpu = false
T = 8.0 # simulation duration in [periods]
nx = 20
dx = 1.0 / nx # pixel resolution in [wavelengths]

l = 2 # simulation domain lxlxl box
sz0 = nx .* (l, l, l)

boundaries = [PEC(-3)] # ground plane on -z, unspecified boundaries default to PML
monitors = [
    # (center, dimensions, normal)
    Monitor([l / 2, l / 2, 1], [1, 1, 0]; normal=[0, 0, 1]),
]
sources = [
    # (signal, center, dimensions)
    Source(t -> cos(2π * t), [l / 2, l / 2, 0.125], [0, 0, 0.25]; Jz=1),
]

configs = maxwell_setup(boundaries, sources, monitors, dx, sz0; F,)
@unpack μ, σ, σm, ϵ, dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

p = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_staggering, p)

# move to gpu
if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, p, field_padding, source_instances = gpu.((u0, p, field_padding, source_instances))
end

# run simulation
@showtime u = accumulate(0:dt:T, init=u0) do u, t
    maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances)
end
y = [power_flux.(u, (m,),) for m = monitor_instances]

# move back to cpu for plotting
if dogpu
    u, p, field_padding, source_instances = cpu.((u, p, field_padding, source_instances))
end

# make movie, 
Ez = get.(u, :Ez)
ϵEz = get(p, :ϵEz)
dir = @__DIR__
recordsim("$dir/$(name).mp4", Ez, y;
    dt,
    field=:Ez,
    monitor_instances,
    source_instances,
    geometry=ϵEz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name Ez"),
    axis2=(; title="monitor powers"),
)
