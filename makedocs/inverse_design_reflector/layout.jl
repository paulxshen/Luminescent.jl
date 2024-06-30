using BSON: @save, @load
using UnPack, Luminescent, GLMakie
include("../../scripts/modes.jl")

F = Float32
dx = 0.05
λ = 1.55f0
ϵ1 = ϵbase = ϵclad = 2.25
ϵ2 = ϵcore = 12.25

wwg = 0.4
hwg = 0.2
lwg = 1.0
ld = 2
hbase = hclad = wm = hm = lm = 0.1
d = 0.25
wd = 2wwg + d
l = lwg + ld + lm
w = wd + 2wm
h = hwg + 2hm

tol = 1.0f-4
function (ε::εtype)(x, y,)
    # lm=()
    if (wm - tol < x < wwg + wm - tol) && (hm - tol < y < hwg + hm - tol)
        return (ϵ2, 0.0, 0.0, ϵ2, ϵ2)
    end

    return (ϵ1, 0.0, 0.0, ϵ1, ϵ1)
end
modes = rectmodes(dx, λ, wwg, hwg, ϵ1, ϵ2, wm, hm,)

function (ε::εtype)(x, y,)
    if (wm - tol < x < wwg + wm - tol) && (hm - tol < y < 4wwg + hm - tol)
        return (ϵ2, 0.0, 0.0, ϵ2, ϵ2)
    end

    return (ϵ1, 0.0, 0.0, ϵ1, ϵ1)
end
sig = linemodes(dx, λ, wwg, ϵ1, ϵ2, wm,)

y = wm + wwg / 2
ports = [
    (; c=[lm, y], lb=[0, -wwg / 2 - wm], ub=[0, wwg / 2 + wm], n=[1, 0]),
    (; c=[lwg - 0.4, d + wm + 1.5wwg], lb=[0, -wwg / 2 - wm], ub=[0, wwg / 2 + wm], n=[-1, 0]),
]
sources = [
    merge((; c=[0, y], n=[1, 0]), sig)
]
designs = [
    (; o=[lwg, wm], L=[ld, wd,],)
]

wwg_, hwg_, lwg_, ld_, wd_, l_, w_, h_, lm_, wm_, d_ =
    round.(Int, [wwg, hwg, lwg, ld, wd, l, w, h, lm, wm, d] ./ dx)

# static_mask = zeros(Int, l_ .+ 1, w_ .+ 1)
static_mask = zeros(Int, l_, w_)
wg = ones(lwg_, wwg_)
place!(static_mask, wg, [1, wm_ + 1],)
place!(static_mask, wg, [1, d_ + wm_ + wwg_ + 1],)
# static_mask[1:lwg_.+1, (w_-wm_-wd_÷2)-wwg_÷2+1:(w_-wm_-wd_÷2)+wwg_÷2+1] .= 1
# static_mask[(l_-lm_-ld_÷2)-wwg_÷2+1:(l_-lm_-ld_÷2)+wwg_÷2+1, 1:lwg_.+1,] .= 1
heatmap(static_mask) |> display
@save "$(@__DIR__)/layout.bson" static_mask sources ports designs dx λ ϵbase ϵclad ϵcore hbase hwg hclad modes l w