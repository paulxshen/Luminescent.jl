using UnPack, LinearAlgebra, Statistics, Random, Jello, DataStructures, FileIO, Porcupine, Dates, NPZ, DataStructures, JSON, Images, BSON, Flux, Zygote, CairoMakie, ArrayPadding, Permutations
using Porcupine: keys, values, fmap, first, ⊙, trim, round, floor, ceil, invperm, permutedims, dict
using Flux: mae, Adam, @functor
using Zygote: withgradient, Buffer, ignore_derivatives, @ignore_derivatives
# using BSON: @save, @load, load
include("core/utils.jl")
include("core/constants.jl")
include("core/modes.jl")
include("core/update.jl")
include("core/boundaries.jl")
include("core/sources.jl")
include("core/monitors.jl")
include("core/geometry.jl")
include("core/gpu.jl")

include("sim/setup.jl")
include("sim/solve.jl")

# include("photonics.jl")
include("format.jl")
# include("dispersion.jl")
# include("c.jl")
include("snapshot.jl")

include("pic/utils.jl")
include("pic/run.jl")
# using AbbreviatedStackTraces

# include("main.jl")