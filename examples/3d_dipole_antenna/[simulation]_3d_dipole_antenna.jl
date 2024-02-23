"""
simulation of quarter wavelength antenna above ground plane
"""

using UnPack, LinearAlgebra, GLMakie
# using FDTDEngine,FDTDToolkit
dir = pwd()
include("$(dir)/src/main.jl")
include("$dir/../FDTDToolkit.jl/src/main.jl")

dogpu = true
# dogpu = false

name = "3d_quarter_wavelength_antenna"
T = 8.0f0 # simulation duration in [periods]
nx = 16
dx = 1.0f0 / nx # pixel resolution in [wavelengths]

l = 2
sz0 = nx .* (l, l, l)

boundaries = [PEC(-3)] # ground plane on -z, unspecified boundaries default to PML
n =
    monitors = [
        # Monitor([.5+l / 2, 0, 1], [0,1, 1,], [1,0,0]),
        Monitor([l / 2, l / 2, 1], [1, 1, 0], [0, 0, 1]),
    ]
sources = [
    Source(t -> cos(F(2π) * t), [l / 2, l / 2, 0.125f0], [0, 0, 0.25f0]; Jz=1),
]

configs = setup(boundaries, sources, monitors, dx, sz0; T)
@unpack μ, σ, σm, ϵ, dt, geometry_padding, geometry_splits, field_padding, source_instances, monitor_instances, u0, = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)


# run simulation
t = 0:dt:T
u = [[similar.(a) for a = u0] for t = t]
u[1] = u0
if dogpu# &&
    using CUDA, Flux
    @assert CUDA.functional()
    u, p, t, field_padding, source_instances = gpu.((u, p, t, field_padding, source_instances))
end

@showtime reduce(
    (u, (u1, t)) -> step!(u1, u, p, t, dx, dt, field_padding, source_instances),
    zip(u[2:end], t[2:end]),
    init=u[1])
y = [power.((m,), u) for m = monitor_instances]

if dogpu# &&
    u, p, t, field_padding, source_instances = cpu.((u, p, t, field_padding, source_instances))
end
Ez = map(u) do u
    u[1][3]
end
dir = @__DIR__
recordsim("$dir/$(name)_nres_$nx.mp4", collect.(Ez), collect.(y);
    dt,
    monitor_instances,
    source_instances,
    geometry=p[1][3],
    elevation=30°,
    axis1=(; title="$name\nEz"),
    axis2=(; title="monitor powers"),
)
