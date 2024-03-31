#=
We simulate a quarter wavelength antenna above conductor ground plane and compute its nearfield radiation pattern
=#

using UnPack, LinearAlgebra, GLMakie, CoordinateTransformations
using GLMakie: volume
using Luminescent, LuminescentVisualization

# if running directly without module # hide
# include("$(pwd())/src/main.jl") # hide
# include("$(pwd())/../LuminescentVisualization.jl/src/main.jl") # hide

name = "quarter_wavelength_antenna"
F = Float32
dogpu = false
T = 8.0 # simulation duration in [periods]
nx = 20
dx = 1.0 / nx # pixel resolution in [wavelengths]

l = 2 # simulation domain lxlxl box
sz = nx .* (l, l, l)
ϵ = ones(F, sz)
μ = ones(F, sz)
σ = zeros(F, sz)
σm = zeros(F, sz)

#=
Set Spherical monitor centered on ground. Portions outside domain eg bottom hemisphere are automatically discarded
=#
boundaries = [PEC(-3)] # ground plane on -z, unspecified boundaries default to PML
monitors = [
    # (center, radius)
    SphereMonitor([l / 2, l / 2, 0], 1),
]
sources = [
    # (signal, center, dimensions)
    Source(t -> cos(2π * t), [l / 2, l / 2, 0.125], [0, 0, 0.25]; Jz=1),
]

configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F,)
@unpack dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

p = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_staggering, p)

# move to gpu
if dogpu
    using Flux
    # using CUDA
    # @assert CUDA.functional()
    u0, p, field_padding, source_instances = gpu.((u0, p, field_padding, source_instances))
end

#=
We run simulation as an `accumulate` loop. `maxwell_update!` applies Maxwells equations as staggered time stepping on E, H. It's mutating so a copy is made in order to save sequence of states
=#
@showtime u = accumulate(0:dt:T, init=u0) do u, t
    maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances)
end

# move back to cpu for plotting
if dogpu
    u, p, field_padding, source_instances = cpu.((u, p, field_padding, source_instances))
end

#=
Plot nearfield Poynting flux thru our Spherical monitor integrated for 1 period
=#
nt = round(Int, 1 / dt)
r = dt * sum(flux.(u[end-nt+1:end], (monitor_instances[1],),))

_, θ, ϕ = eachrow(sphcoords(monitors[1])[:, inbounds(monitor_instances[1])])
cfs = CartesianFromSpherical()
rvecs = cfs.(splat(Spherical).(zip(r, F.(θ), F.(ϕ))))

fig = Figure()
ax = Axis3(fig[1, 1])
plot!(ax, getindex.(rvecs, 1), getindex.(rvecs, 2), getindex.(rvecs, 3),)
display(fig)
save("antennapattern.png", fig)

#=
![](assets/antennapattern.png)
=#

#=
Ready, set, action! We make movie, 
=#
Ez = field.(u, :Ez)
ϵEz = field(p, :ϵEz)
dir = @__DIR__
recordsim("$dir/$(name).mp4", Ez, ;
    dt,
    field=:Ez,
    monitor_instances,
    source_instances,
    geometry=ϵEz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name Ez"),
    # axis2=(; title="monitor powers"),
)

#=
![](assets/quarter_wavelength_antenna.mp4)
=#