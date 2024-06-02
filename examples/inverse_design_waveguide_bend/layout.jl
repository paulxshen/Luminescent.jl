using BSON: @save, @load
include("../../scripts/modes.jl")

F = Float32
dx = 0.05
λ = 1.55f0
ϵ1 = ϵbase = ϵclad = 2.25
ϵ2 = ϵcore = 12.25

wwg = 0.4
hwg = 0.2
lwg = 1.5
ld = 1.6
wd = 1.6
hbase = hclad = wm = hm = lm = 0.25
l = lwg + ld + lm
w = lwg + wd + wm
h = hwg + 2hm

# wwg = 0.5
# hwg = 0.25


sig = rectmodes("$(@__DIR__)/modes.bson", dx, λ, wwg, hwg, wm, hm, ϵ1, ϵ2)
ports = [
    (; c=[lwg - lm, w - wm - wd / 2], lb=[0, -wwg / 2 - wm], ub=[0, wwg / 2 + wm], n=[1, 0]),
    (; c=[l - lm - ld / 2, wm], lb=[-wwg / 2 - wm, 0], ub=[wwg / 2 + wm, 0], n=[0, -1]),
]
sources = [
    merge((; c=[0, w - wm - wd / 2], n=[1, 0]), sig)
]
designs = [
    (; o=[lwg, lwg], L=[ld, wd,],)
]

wwg_, hwg_, lwg_, ld_, wd_, l_, w_, h_, lm_, wm_ =
    round.(Int, [wwg, hwg, lwg, ld, wd, l, w, h, lm, wm] ./ dx)

static_mask = zeros(Int, l_ .+ 1, w_ .+ 1)
wg = ones(lwg_, wwg_)
place!(static_mask, wg, [1, (w_ - wm_ - wd_ ÷ 2) - wwg_ ÷ 2 + 1],)
place!(static_mask, wg', [(l_ - lm_ - ld_ ÷ 2) - wwg_ ÷ 2 + 1, 1],)
# static_mask[1:lwg_.+1, (w_-wm_-wd_÷2)-wwg_÷2+1:(w_-wm_-wd_÷2)+wwg_÷2+1] .= 1
# static_mask[(l_-lm_-ld_÷2)-wwg_÷2+1:(l_-lm_-ld_÷2)+wwg_÷2+1, 1:lwg_.+1,] .= 1

@save "$(@__DIR__)/layout.bson" static_mask sources ports designs dx λ ϵbase ϵclad ϵcore hbase hwg hclad