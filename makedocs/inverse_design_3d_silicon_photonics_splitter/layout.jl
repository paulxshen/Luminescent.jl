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
o = [lwg, wm] / λ
v = F.(v) ./ λ
# L = [l, w, h]
ports = [[0.2f0, [(w - wwg) / 2, (w + wwg) / 2]], [l, [w - wm - wwg, w - wm]], [l, [wm, wm + wwg]]] / λ
sources = [[0, w / 2]] / λ
designs = [[lwg:lwg+ld, wm:wm+ld]] / λ


dx = 1.0f0 / nres
wwg, hwg, lwg, ld, hclad, hbox, l, w, h, wm = round.(Int, v ./ dx)
base = zeros(Int, l, w)
base[1:lwg, (w-wwg)÷2+1:(w+wwg)÷2] .= 1
base[end-lwg+1:end, wm+1:wm+wwg] .= 1
base[end-lwg+1:end, end-wm-wwg+1:end-wm] .= 1
# design = zeros(F, l, w)
# design[lwg:lwg+ld, wwg:wwg+ld] = 1

layout = (; Pair.(k, v)..., base, sources, ports, designs, dx, o, nres)
# @save "$dir/layout.bson" layout