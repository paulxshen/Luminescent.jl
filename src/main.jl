using UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures
using Porcupine
using ChainRules: ignore_derivatives
using Jello, ArrayPadding
using Zygote: bufferfrom, Buffer
using Zygote
F = Float32
Random.seed!(1)
include("utils.jl")
# include("../../ArrayPadding.jl/src/main.jl")
# include("../../Porcupine.jl/src/del.jl")
include("maxwell_update.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("maxwell_setup.jl")
include("geometry.jl")

# ]add UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures,Porcupine,Jello, ArrayPadding,Zygote, GPUArraysCore