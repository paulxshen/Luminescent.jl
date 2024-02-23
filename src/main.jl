using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Functors, DataStructures
using Porcupine
using Jello, ArrayPadding
using Zygote: bufferfrom, Buffer
using Zygote
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

