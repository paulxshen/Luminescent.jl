using BSON: @save, @load
using UnPack, Luminescent
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
wd = wwg
# wd = 2
hsub = hclad = wm = hm = lm = 0.25
l = 2lwg + ld
w = wd + 3wm
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

ports = [
    (; c=[lm, 2wm + wwg / 2], lb=[0, -wwg / 2 - wm], ub=[0, wwg / 2 + wm], n=[1, 0]),
    (; c=[l - lm, 2wm + wwg / 2], lb=[0, -wwg / 2 - wm], ub=[0, wwg / 2 + wm], n=[1, 0]),
]
signals = [
    merge((; c=[0, 2wm + wwg / 2], n=[1, 0]), sig)
]
designs = [
    (; o=[lwg, 2wm], L=[ld, wd,],)
]

wwg_, hwg_, lwg_, ld_, wd_, l_, w_, h_, lm_, wm_ =
    round.(Int, [wwg, hwg, lwg, ld, wd, l, w, h, lm, wm] ./ dx)

# static_mask = zeros(Int, l_ .+ 1, w_ .+ 1)
static_mask = zeros(Int, l_, w_)
wg = ones(lwg_, wwg_)
# y = (w_ ÷ 2) - wwg_ ÷ 2 + 1
y = 2wm_ + 1
place!(static_mask, wg, [1, y],)
place!(static_mask, wg, [lwg_ + ld_ + 1, y],)
# static_mask[1:lwg_.+1, (w_-wm_-wd_÷2)-wwg_÷2+1:(w_-wm_-wd_÷2)+wwg_÷2+1] .= 1
# static_mask[(l_-lm_-ld_÷2)-wwg_÷2+1:(l_-lm_-ld_÷2)+wwg_÷2+1, 1:lwg_.+1,] .= 1

@save "$(@__DIR__)/layout.bson" static_mask signals ports designs dx λ ϵbase ϵclad ϵcore hsub hwg hclad modes