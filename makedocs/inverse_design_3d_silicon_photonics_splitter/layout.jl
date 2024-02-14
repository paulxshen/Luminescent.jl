using BSON: @save, @load
F = Float32
dir = "examples/utils"

wwg = 0.5
hwg = 0.22
lwg = 1
ld = 2
hclad = hbox = 0.2
wm = 0.5wwg
l = 2lwg + ld
w = 2wm + ld
h = hwg + hbox + hclad
layout = (; wwg, hwg, lwg, ld, hclad, hbox, l, w, h, wm)
v = values(layout)
k = keys(layout)

λ = 1.55f0
c = [lwg, wm] / λ
v = F.(v) ./ λ
# L = [l, w, h]
ports = [
    (c = [0.2f0, w / 2] / λ),
    (c = [l, w - wm] / λ),
    (c = [l, wm] / λ),
]
sources = [
    (c = [0, w / 2] / λ)
]
designs = [
    (c=[lwg, wm] / λ, L=[ld, ld,] / λ)
]


dx = 1.0f0 / nres
wwg, hwg, lwg, ld, hclad, hbox, l, w, h, wm = round.(Int, v ./ dx)
base = zeros(Int, l, w)
base[1:lwg, (w-wwg)÷2+1:(w+wwg)÷2] .= 1
base[end-lwg+1:end, wm+1:wm+wwg] .= 1
base[end-lwg+1:end, end-wm-wwg+1:end-wm] .= 1
# design = zeros(F, l, w)
# design[lwg:lwg+ld, wwg:wwg+ld] = 1

layout = (; Pair.(k, v)..., base, sources, ports, designs, dx, nres)
@save "$(@__DIR__)/layout.bson" layout