using UnPack, LinearAlgebra, Statistics, Random, Jello, DataStructures, FileIO, Porcupine, Dates, NPZ, DataStructures, JSON, Flux, Zygote, CairoMakie, ArrayPadding, Permutations, Functors, Optimisers, ChainRulesCore, SparseArrays, Interpolations
using Meshes, GeoIO
using Meshes: Point, boundingbox
# using GeometryBasics, MeshIO
using Porcupine: keys, values, pairs, fmap, ⊙, trim, round, floor, ceil, invperm, permutedims, dict, cpu, gpu
using Flux: mae, Adam
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
include("core/mesh.jl")

include("sim/setup.jl")
include("sim/solve.jl")

include("format.jl")
# include("dispersion.jl")
# include("c.jl")
include("snapshot.jl")

include("pic/utils.jl")
include("pic/run.jl")

include("gen/run.jl")

include("ops.jl")
include("del.jl")
# include("main.jl")

include("sparams.jl")
include("plots.jl")
include("gpu.jl")

function julia_main()::Cint
    if !isempty(ARGS)
        picrun(ARGS[1])
    end
    return 0
end