using Porcupine: keys, values, fmap, first, ⊙, trim
using Humanize: digitsep
using UnPack, LinearAlgebra, Statistics, Random, Jello, Functors, DataStructures, GPUArraysCore, TrackedFloats
using Zygote
using CairoMakie

using Porcupine
using ArrayPadding
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
# include("c.jl")
include("snapshot.jl")

using Dates, DataStructures, JSON, Images, BSON, Flux, CUDA, GPUArraysCore
using Flux: mae, Adam
using Zygote: withgradient, Buffer, ignore_derivatives
using BSON: @save, @load, load
include("gf.jl")
# using AbbreviatedStackTraces

# include("main.jl")