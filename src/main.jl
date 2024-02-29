using UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Functors, DataStructures
using Porcupine
using ChainRules: ignore_derivatives
using Jello, ArrayPadding
using Zygote: bufferfrom, Buffer
using Zygote
using GPUArraysCore
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

