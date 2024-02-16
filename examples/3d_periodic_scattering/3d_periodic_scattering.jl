"""
simulation of plane wave scattering on periodic array of dielectric spheres
"""

using UnPack, LinearAlgebra, GLMakie
# using FDTDEngine,FDTDToolkit
dir = pwd()
include("$(dir)/src/main.jl")
include("$dir/../FDTDToolkit.jl/src/main.jl")


name = "3d_scattering"
T = 8.0f0 # simulation duration in [periods]
nx = 16
dx = 1.0f0 / nx # pixel resolution in [wavelengths]

"geometry"
l = 2 # domain physical size length
sz = nx .* (l, l, l) # domain voxel dimensions
ϵ1 = ϵmin = 1 #
ϵ2 = 2.25f0 # 
b = F.([norm(v .- sz ./ 2) < 0.5 / dx for v = Base.product(Base.oneto.(sz)...)]) # sphere
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

"setup"
boundaries = [Periodic(2), Periodic(3)]# unspecified boundaries default to PML
sources = [
    PlaneWave(t -> cos(F(2π) * t), -1; Jz=1) # Jz excited plane wave from -x plane (eg -1)
]
n = [1, 0, 0] # normal 
δ = 0.2f0 # margin
# A = (l - δ)^2
lm = 1 # monitor side length
monitors = [
    Monitor([δ, l / 2, l / 2], [0, lm, lm], n,), # (center, dimensions, normal)
    Monitor([l - δ, l / 2, l / 2], [0, lm, lm], n,),
]
configs = setup(boundaries, sources, monitors, dx, sz; ϵmin, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_instances, monitor_instances, u0, = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)


# run simulation
t = 0:dt:T
u = [[similar.(a) for a = u0] for t = t]
u[1] = u0
@showtime reduce(
    (u, (u1, t)) -> step!(u1, u, p, t, dx, dt, field_padding, source_instances),
    zip(u[2:end], t[1:end-1]),
    init=u0)
y = hcat([power.((m,), u) for m = monitor_instances]...)

# make movie
Ez = map(u) do u
    u[1][3]
end
ϵz = p[1][3]
dir = @__DIR__

recordsim("$dir/$(name)_nres_$nx.mp4", Ez, y;
    dt,
    monitor_instances,
    source_instances,
    geometry=ϵz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name\nEz"),
    axis2=(; title="monitor powers"),
)
