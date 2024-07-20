using Humanize: digitsep
using UnPack, LinearAlgebra, Random, Jello, StatsBase, ImageTransformations, Meshes, CoordinateTransformations, Functors, DataStructures, GPUArraysCore
# using Meshes: Sphere
using ArrayPadding: left, right
using Zygote: bufferfrom, Buffer, ignore_derivatives
using Zygote, Flux#,NNlib,
using CairoMakie
# Random.seed!(1)

using Porcupine: keys, values, fmap
using Porcupine
using ArrayPadding
# include("../../Porcupine.jl/src/main.jl")
# include("../../ArrayPadding.jl/src/main.jl")
# include("../../Jello.jl/src/main.jl")
include("constants.jl")
include("utils.jl")
include("modes.jl")
include("update.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("geometry.jl")
include("setup.jl")
include("solve.jl")
include("photonics.jl")
# include("gpu.jl")
include("format.jl")
include("dispersion.jl")
include("snapshot.jl")
# ]add UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures,Porcupine,Jello, ArrayPadding,Zygote, GPUArraysCore