# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
Random.seed!(1234)
using CUDA
# prob, sol = genrun("phantom"; Ttrans=3)
# prob, sol = genrun("ant"; Ttrans=10)
using GLMakie
# S, sol, prob = genrun(joinpath("genruns", "ant"), cu)
S, sol, prob = genrun(joinpath("genruns", "wg"), cu)
# S, sol, prob = genrun(joinpath("G:", "My Drive", "reflex", "wg"), cu)
# S, sol, prob = genrun(joinpath(raw"G:\My Drive\reflex\wg"), cu)
# vis(sol, prob, :Ey)

# a = sol.u(:Ey)
# volume(abs.(a))
# extrema(a)
# abs(sum(sol.um[1][1].Ey) / sum(sol.um[1][1].Hx)) |> println

# using GLMakie
# a = 10prob.source_instances[1].sigmodes[1][2](2)
# p = prob.geometry
# p = pad_geometry(p, prob.geometry_padvals, prob.geometry_padamts)
# b = p.ϵ
# b /= maximum(b)
# volume(a + b)

# volume(sources[1].mask)
# prob.source_instances[1].sigmodes[1][1].([0, 0.25, 0.5, 0.75])

# a = _p.ϵ(1) /
#     # a = _p.σ(1)
#     extrema(a) |> println
# volume(a)
# heatmap(a[:, :, round(Int, end ÷ 2)])

# a = sol.u(:Ey) |> cpu
# volume(abs.(a))
# v = maximum(abs, a)
# a /= v
# # a = abs.(a) / v
# # a = min.(a, 0.01)
# @show extrema(a)
# r = 0.1
# colorrange = (-r, r)
# # colormap = cgrad([RGBAf(x, 0, 1 - x, (2abs(x - 0.5))^0.3) for x = range(0, 1, length=11)])
# colormap = cgrad([:blue, RGBAf(0, 0, 0, 0), :red])
# volume(a; colormap, colorrange, algorithm=:absorption)

# # volume(sources[1].mask)
# # volume(_aa)
0