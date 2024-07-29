using Humanize: digitsep
using UnPack, LinearAlgebra, Statistics, Random, Jello, Functors, DataStructures, GPUArraysCore
# using Meshes: Sphere
using ArrayPadding: left, right
using Zygote
using CairoMakie
# Random.seed!(1)ImageIO,ImageTransformations,

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
include("gpu.jl")
include("format.jl")
include("dispersion.jl")
include("snapshot.jl")

using Dates, DataStructures, JSON, Images, BSON, Flux, CUDA, GPUArraysCore
using Flux: mae, Adam
using Zygote: withgradient, Buffer, ignore_derivatives
using BSON: @save, @load, load
# using AbbreviatedStackTraces

include("gf.jl")
# include("main.jl")