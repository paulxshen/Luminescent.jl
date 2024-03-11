using BSON: @save, @load
F = Float32
dir = "examples/utils"

nx = 31
wwg = 0.5
hwg = 0.25
lwg = 0.7
lwg1 = 0
ldx = 1.5
ldy = 2
hsub = hclad = lm = wm = 0.25
l = 2lwg + lwg1 + ldx
w = 2wm + ldy
h = hwg + hsub + hclad

λ = 1.55f0

c = [lwg, wm]
n = [1, 0]
δ = 0.1
ports = [
    (; c=[lwg, w / 2], n),
    (; c=[l - δ, w - wm - wwg], n),
    (; c=[l - δ, wm + wwg], n),
]
signals = [
    (; c=[0, w / 2])
]
designs = [
    (; o=[lwg + lwg1, wm], L=[ldx, ldy,], c=[lwg + lwg1 + ldx / 2, wm + ldy / 2])
]


dx = 0.05

wwg, hwg, lwg, ldx, ldy, hclad, hsub, l, w, h, wm =
    round.(Int, [wwg, hwg, lwg, ldx, ldy, hclad, hsub, l, w, h, wm] ./ dx)

base = zeros(Int, l .+ 1, w .+ 1)
base[1:lwg.+1, (w-wwg)÷2+1:(w+wwg)÷2+1] .= 1
base[end-lwg:end, wm+1:wm+wwg.+1] .= 1
base[end-lwg:end, end-wm-wwg:end-wm] .= 1

@save "$(@__DIR__)/layout.bson" base signals ports designs dx nx