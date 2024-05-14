using UnPack, LinearAlgebra, Random, Jello, StatsBase, ImageTransformations, Meshes, CoordinateTransformations, Functors, DataStructures, ArrayPadding
using Meshes: Sphere
using ChainRules: ignore_derivatives
using ArrayPadding: left, right
using Zygote: bufferfrom, Buffer
using Zygote, Flux
# Random.seed!(1)

# using Porcupine: keys, values
# using  Porcupine
include("../../Porcupine.jl/src/main.jl")
# include("../../Jello.jl/src/main.jl")

include("constants.jl")
include("utils.jl")
include("modes.jl")
include("maxwell_update.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("geometry.jl")
include("maxwell_setup.jl")
include("photonics.jl")
# ]add UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures,Porcupine,Jello, ArrayPadding,Zygote, GPUArraysCore