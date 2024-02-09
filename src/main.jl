using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Functors, DataStructures
using Jello, ArrayPadding, Porcupine
using Zygote: bufferfrom, Buffer
F = Float32
Random.seed!(1)
include("utils.jl")
# include("../../ArrayPadding.jl/src/main.jl")
# include("../../Porcupine.jl/src/del.jl")
include("maxwell.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("fdtd.jl")
include("geometry.jl")

