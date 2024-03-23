using BSON: @save, @load
include("../../scripts/modes.jl")

F = Float32
dx = 0.05
λ = 1.55f0
ϵ1 = ϵsub = ϵclad = 2.25
ϵ2 = ϵcore = 12.25

wwg = 0.4
hwg = 0.2
lwg = 1.5
ld = 1.6
wd = 1.6
hsub = hclad = wm = hm = lm = 0.25
l = lwg + ld + lm
w = lwg + wd + wm
h = hwg + 2hm

# wwg = 0.5
# hwg = 0.25


sig = savemodes("$(@__DIR__)/modes.bson", dx, λ, wwg, hwg, wm, hm, ϵ1, ϵ2)
ports = [
    (; c=[lwg - lm, w - wm - wd / 2], lb=[0, -wwg / 2 - δ], ub=[0, wwg / 2 + δ], n=[1, 0]),
    (; c=[l - lm - ld / 2, wm], lb=[-wwg / 2 - δ, 0], ub=[wwg / 2 + δ, 0], n=[0, -1]),
]
signals = [
    merge((; c=[0, w - wm - wd / 2], n=[1, 0]), sig)
]
designs = [
    (; o=[lwg, lwg], L=[ld, wd,],)
]

wwg_, hwg_, lwg_, ld_, wd_, l_, w_, h_, lm_, wm_ =
    round.(Int, [wwg, hwg, lwg, ld, wd, l, w, h, lm, wm] ./ dx)

static_mask = zeros(Int, l_ .+ 1, w_ .+ 1)
static_mask[1:lwg_.+1, (w_-wm_-wd_÷2)-wwg_÷2+1:(w_-wm_-wd_÷2)+wwg_÷2+1] .= 1
static_mask[(l-lm_-ld_÷2)-wwg_÷2+1:(l_-lm_-ld_÷2)+wwg_÷2+1, 1:lwg_.+1,] .= 1

@save "$(@__DIR__)/layout.bson" static_mask signals ports designs dx λ ϵsub ϵclad ϵcore hsub hwg hclad