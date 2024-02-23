using BSON: @save, @load
F = Float32
dir = "examples/utils"

wwg = 0.5
hwg = 0.22
lwg = 1
ld = 1.5
hclad = hbox = 0.2
wm = 0.5wwg
l = 2lwg + ld
w = 2wm + ld
h = hwg + hbox + hclad

λ = 1.55f0
layout = (; wwg, hwg, lwg, ld, hclad, hbox, l, w, h, wm)
v = collect(values(layout)) / λ
k = keys(layout)

c = [lwg, wm]
v = F.(v)
n = [1, 0]
ports = [
    (; c=[0.2f0, w / 2], n),
    (; c=[l, w - wm - wwg / 2], n),
    (; c=[l, wm + wwg / 2], n),
]
sources = [
    (; c=[0, w / 2])
]
designs = [
    (; o=[lwg, wm], L=[ld, ld,], c=[lwg + ld / 2, wm + ld / 2])
]


dx = 1.0f0 / nx
wwg, hwg, lwg, ld, hclad, hbox, l, w, h, wm = round.(Int, v ./ dx)

base = zeros(Int, l .+ 1, w .+ 1)
base[1:lwg.+1, (w-wwg)÷2+1:(w+wwg)÷2+1] .= 1
base[end-lwg:end, wm+1:wm+wwg.+1] .= 1
base[end-lwg:end, end-wm-wwg:end-wm] .= 1

layout = (; Pair.(k, v)..., base, sources, ports, designs, dx, nx)
@save "$(@__DIR__)/layout.bson" layout