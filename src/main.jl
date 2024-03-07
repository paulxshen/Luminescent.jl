using UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures
using Porcupine
using ChainRules: ignore_derivatives
using Jello, ArrayPadding
using Zygote: bufferfrom, Buffer
using Zygote
# Random.seed!(1)
include("utils.jl")
# include("../../ArrayPadding.jl/src/main.jl")
# include("../../Porcupine.jl/src/del.jl")
include("maxwell_update.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("geometry.jl")
include("maxwell_setup.jl")

# ]add UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures,Porcupine,Jello, ArrayPadding,Zygote, GPUArraysCore