using Porcupine: keys, values, fmap, first, ⊙, trim, round, floor, ceil
using Humanize: digitsep
using UnPack, LinearAlgebra, Statistics, Random, Jello, Functors, DataStructures, GPUArraysCore, FileIO
using Zygote
using CairoMakie

using Porcupine
using Dates, DataStructures, JSON, Images, BSON, Flux, CUDA, GPUArraysCore
using Flux: mae, Adam
using Zygote: withgradient, Buffer, ignore_derivatives
using BSON: @save, @load, load
using ArrayPadding
include("constants.jl")
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

include("pic/utils.jl")
include("pic/run.jl")
# using AbbreviatedStackTraces

# include("main.jl")