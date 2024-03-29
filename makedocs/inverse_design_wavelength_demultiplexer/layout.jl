using BSON: @save, @load
F = Float32
dir = "examples/utils"

nx = 31
wwg = 0.5
hwg = 0.25
lwg = 0.5
ld = 2
hsub = hclad = lm = wm = 0.25
l = 2lwg + ld
w = 2wm + ld
h = hwg + hsub + hclad

λ = 1.55f0

c = [lwg, wm]
n = [1, 0]
δ = 0.1
ports = [
    (; c=[δ, w / 2], n),
    (; c=[l - δ, w - wm - wwg], n),
    (; c=[l - δ, wm + wwg], n),
]
signals = [
    (; c=[0, w / 2])
]
designs = [
    (; o=[lwg, wm], L=[ld, ld,], c=[lwg + ld / 2, wm + ld / 2])
]


dx = 0.05

wwg, hwg, lwg, ld, hclad, hsub, l, w, h, wm =
    round.(Int, [wwg, hwg, lwg, ld, hclad, hsub, l, w, h, wm] ./ dx)

mask = zeros(Int, l .+ 1, w .+ 1)
mask[1:lwg.+1, (w-wwg)÷2+1:(w+wwg)÷2+1] .= 1
mask[end-lwg:end, wm+1:wm+wwg.+1] .= 1
mask[end-lwg:end, end-wm-wwg:end-wm] .= 1

@save "$(@__DIR__)/layout.bson" mask signals ports designs dx nx