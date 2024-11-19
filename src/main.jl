using UnPack, LinearAlgebra, Statistics, Random, Jello, DataStructures, FileIO, Porcupine, Dates, DataStructures, JSON, Images, BSON, Flux, CUDA, Zygote, CairoMakie, ArrayPadding
using Porcupine: keys, values, fmap, first, ⊙, trim, round, floor, ceil
using Flux: mae, Adam, @functor
using Zygote: withgradient, Buffer, ignore_derivatives
using BSON: @save, @load, load
include("utils.jl")
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