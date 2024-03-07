"""
simulation of  coupling into dielectric slab waveguide using modal source 
"""

using UnPack, LinearAlgebra, GLMakie
using BSON: @load
using Luminescent, LuminescentVisualization

# dir = pwd()
# include("$(dir)/src/main.jl")
# include("$(dir)/scripts/startup.jl")
# include("$dir/../LuminescentVisualization.jl/src/main.jl")

dogpu = true
name = "slab_waveguide"
T = 14 # simulation duration in [periods]

# load mode profile and waveguide dimensions from results of external mode solver 
@load "$(@__DIR__)/modes.bson" modes lb ub λ dx hsub wwg hwg hclad w h ϵsub ϵclad ϵwg
ϵmin = ϵclad
hsub, wwg, hwg, hclad, w, dx, ub, lb = [hsub, wwg, hwg, hclad, w, dx, ub, lb] / λ

# geometry
l = 2
mask = zeros(round.(Int, (l, w) ./ dx))

# waveguide slab
mask[:, round(Int, (w / 2 - wwg / 2) / dx)+1:round(Int, (w / 2 + wwg / 2) / dx)+1] .= 1
ϵ = sandwich(mask, round.(Int, [hsub, hwg, hclad] / dx), [ϵsub, ϵwg, ϵclad])
sz = size(ϵ)


# modal source
@unpack Ex, Ey, Ez, = modes[1]
Jy, Jz, Jx = map([Ex, Ey, Ez]) do a
    reshape(a, 1, size(a)...)
end
# GLMakie.volume(real(Jy))
c = [0, w / 2, hsub]
lb = [0, lb...]
ub = [0, ub...]
sources = [Source(t -> cispi(2t), c, lb, ub; Jx, Jy, Jz)]

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
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; ϵmin,)
@unpack μ, σ, σm, dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_staggering; ϵ, μ, σ, σm)

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
v = [power_flux.((m,), u) for m = monitor_instances]

# move back to cpu for plotting
if dogpu
    u, p, field_padding, source_instances = cpu.((u, p, field_padding, source_instances))
end

# make movie, 
i = 2
Ey = field.(u, :Ey)
ϵEy = p[1][i]
dir = @__DIR__
recordsim("$dir/$(name).mp4", Ey, v;
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
