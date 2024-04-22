"""
simulation of  coupling into dielectric slab waveguide using modal source . We do a 2d TE simulation followed by 3d. The former is a fast approximation of the later
"""

using UnPack, LinearAlgebra, GLMakie
using BSON: @load
# using Luminescent, LuminescentVisualization

dir = pwd()
include("$(dir)/src/main.jl")
include("$(dir)/scripts/startup.jl")
include("$dir/../LuminescentVisualization.jl/src/main.jl")

name = "slab_waveguide"
dogpu = true
F = Float32

# load mode profile and waveguide dimensions from results of external mode solver 
@load "$(@__DIR__)/modes.bson" modes lb ub λ dx hsub wwg hwg hclad w h ϵbase ϵclad ϵcore
ϵmin = ϵclad
hsub, wwg, hwg, hclad, w, dx, ub, lb = [hsub, wwg, hwg, hclad, w, dx, ub, lb] / λ
dx = F(dx)

# geometry
l = 2
mask = zeros(round.(Int, (l, w) ./ dx))

# waveguide slab
mask[:, round(Int, (w / 2 - wwg / 2) / dx)+1:round(Int, (w / 2 + wwg / 2) / dx)+1] .= 1
sz = size(mask)

ϵ = mask * ϵcore + (1 .- mask) * ϵclad
μ = 1
σ = zeros(F, sz)
σm = zeros(F, sz)
T = 2 + l * sqrt(ϵcore) # simulation duration in [periods]

# plot(abs.(source_instances[1]._g[:Jx]))
# modal source
@unpack Ex, Ey, Ez, = modes[1]
Jy, Jx = map([Ex, Ez]) do a
    transpose(sum(a, dims=2))
end
Jy, Jx = [Jy, Jx] / maximum(maximum.(abs, [Jy, Jx]))
c = [0, w / 2,]
lb_ = [0, lb[1]]
ub_ = [0, ub[1]]
sources = [Source(t -> cispi(2t), c, lb_, ub_; Jx, Jy,)]

# monitors
normal = [1, 0] # normal 
δ = 0.1 / λ # margin
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor([l / 2, w / 2,], [0, -wwg / 2 - δ,], [0, wwg / 2 + δ,]; normal,),
    Monitor([l - δ, w / 2,], [0, -wwg / 2 - δ,], [0, wwg / 2 + δ,]; normal,),
]

# maxwell_setup
boundaries = []# unspecified boundaries default to PML
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; ϵmin, F)
@unpack dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

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
v = [power.(u, (m,),) for m = monitor_instances]

# move back to cpu for plotting
if dogpu
    u, p, field_padding, source_instances = cpu.((u, p, field_padding, source_instances))
end

# make movie, 
Ey = field.(u, :Ey)
ϵEy = field(p, :ϵEy)
# heatmap(ϵEy)
dir = @__DIR__
recordsim("$dir/2d_$(name).mp4", Ey, v;
    dt,
    field=:Ey,
    monitor_instances,
    source_instances,
    geometry=ϵEy,
    # geometry=mask,
    elevation=30°,
    playback=1,
    axis1=(; title="$name Ey"),
    axis2=(; title="monitor powers"),
)
# error()

# begin 3d simulation
ϵ = sandwich(mask, round.(Int, [hsub, hwg, hclad] / dx), [ϵbase, ϵcore, ϵclad])
sz = size(ϵ)
μ = 1
σ = zeros(F, sz)
σm = zeros(F, sz)

# modal source
@unpack Ex, Ey, Ez, = modes[1]
Jy, Jz, Jx = map([Ex, Ey, Ez]) do a
    reshape(a, 1, size(a)...)
end
Jz, Jy, Jx = [Jz, Jy, Jx] / maximum(maximum.(abs, [Jz, Jy, Jx]))
# GLMakie.volume(real(Jy))
c = [0, w / 2, hsub]
lb_ = [0, lb...]
ub_ = [0, ub...]
sources = [Source(t -> cispi(2t), c, lb_, ub_; Jx, Jy, Jz)]

# monitors
normal = [1, 0, 0] # normal 
δ = 0.1 / λ # margin
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor([l / 2, w / 2, hsub], [0, -wwg / 2 - δ, -δ], [0, wwg / 2 + δ, hwg + δ]; normal,),
    Monitor([l - δ, w / 2, hsub], [0, -wwg / 2 - δ, -δ], [0, wwg / 2 + δ, hwg + δ]; normal,),
]

# maxwell_setup
boundaries = []# unspecified boundaries default to PML
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; ϵmin, F)
@unpack dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

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
v = [power.(u, (m,),) for m = monitor_instances]

# move back to cpu for plotting
if dogpu
    u, p, field_padding, source_instances = cpu.((u, p, field_padding, source_instances))
end

# make movie, 
Ey = field.(u, :Ey)
ϵEy = field.(p, :ϵEy)
dir = @__DIR__
recordsim("$dir/3d_$(name).mp4", Ey, v;
    dt,
    field=:Ey,
    monitor_instances,
    source_instances,
    geometry=ϵEy,
    elevation=30°,
    playback=1,
    axis1=(; title="$name Ey"),
    axis2=(; title="monitor powers"),
)
