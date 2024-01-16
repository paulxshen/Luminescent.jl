using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Flux, Functors, Optim, DataStructures
using Flux: withgradient
using Optim: Options, minimizer
using Zygote: ignore, Buffer
using Jello, ArrayPadding,Porcupine

F = Float32
Random.seed!(1)
include("utils.jl")
include("maxwell.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("fdtd.jl")

