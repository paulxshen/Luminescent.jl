using UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Meshes, CoordinateTransformations, Functors, DataStructures
using Porcupine: keys, values
using Meshes: Sphere
using ChainRules: ignore_derivatives
using ArrayPadding, Porcupine
using ArrayPadding: left, right
using Zygote: bufferfrom, Buffer
using Zygote, Flux
using Jello
# Random.seed!(1)
# include("../../Porcupine.jl/src/main.jl")
# include("../../Jello.jl/src/main.jl")
include("maxwell_update.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("geometry.jl")
include("utils.jl")
include("maxwell_setup.jl")

# ]add UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures,Porcupine,Jello, ArrayPadding,Zygote, GPUArraysCore