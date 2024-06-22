using UnPack, LinearAlgebra, Random, Jello, StatsBase, ImageTransformations, Meshes, CoordinateTransformations, Functors, DataStructures
using Meshes: Sphere
# using ArrayPadding: left, right
using Zygote: bufferfrom, Buffer, ignore_derivatives
using Zygote, Flux#,NNlib
# Random.seed!(1)

# using Porcupine: keys, values
# using  Porcupine
include("../../Porcupine.jl/src/main.jl")
include("../../ArrayPadding.jl/src/main.jl")
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
include("solve.jl")
include("photonics.jl")
# ]add UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures,Porcupine,Jello, ArrayPadding,Zygote, GPUArraysCore