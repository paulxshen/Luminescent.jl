"""
simulation of quarter wavelength antenna above ground plane
"""

using UnPack, LinearAlgebra, GLMakie
# using FDTDEngine,FDTDToolkit
dir = pwd()
include("$(dir)/src/main.jl")
include("$dir/../FDTDToolkit.jl/src/main.jl")

name = "3d_quarter_wavelength_antenna"
T = 8.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]

l = 2
sz0 = nres .* (l, l, l)

boundaries = [PEC(-3)] # unspecified boundaries default to PML
n = [0, 0, 1]
monitors = [
    Monitor([l / 2, l / 2, 0.1f0], [1, 1, 0], n),
    Monitor([l / 2, l / 2, 0.5f0], [1, 1, 0], n),
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
u = similar([u0], length(t))
@showtime u = accumulate!((u, t) -> step!(u, p, t, dx, dt, field_padding, source_instances), u, t, init=deepcopy(u0))
y = hcat([power.((m,), u) for m = monitor_instances]...)

Ez = map(u) do u
    u[1][3]
end
dir = @__DIR__
@showtime recordsim("$dir/$(name)_nres_$nres.mp4", Ez, y;
    dt,
    monitor_instances,
    source_instances,
    geometry=p[1][3],
    elevation=30°,
    axis1=(; title="$name\nEz"),
    axis2=(; title="monitor powers"),
)
