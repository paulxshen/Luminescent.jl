# using CairoMakie, Random
# Random.seed!(1)
# n = 50
# a = rand(Float32, n, n)
# @time heatmap(a)
# @time heatmap(a)
# @time heatmap(a)

using Random, LinearAlgebra
using Jello
l = 100

Random.seed!(1)
contrast = 1
rmin = nothing
m = Blob(l, l; contrast)
