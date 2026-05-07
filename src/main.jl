using UnPack, LinearAlgebra, Statistics, StatsBase, Random, Dates, NPZ, DataStructures, JSON, Zygote, Unitful, NearestNeighbors, Functors, NNlib, GPUArraysCore, GPUArrays, OhMyThreads, Interpolations, ChainRulesCore, SparseArrays, StaticArrays, Meshes, GeoIO, FileIO, Colors, UnPack, Tar
using Meshes: Point, boundingbox, Box
using Zygote: withgradient, Buffer, ignore_derivatives, @ignore_derivatives
# using CUDA
using ArrayPadding, Jello, Porcupine
using Porcupine: rmap, ⊙, dict, cpu, gpu, fmap
using GPUArraysCore: @allowscalar
# using GLMakie
using CairoMakie

# const VERSION = v"2.4.0"
# function Zygote.gradindex(x::Tuple, i::Int)
#     i > length(x) && return NoTangent()
#     x[i]
# end

include("log.jl")
include("core/constants.jl")

include("core/utils.jl")
include("core/update.jl")
include("core/boundaries.jl")
include("core/monitors.jl")
include("core/modes.jl")
include("core/sources.jl")
include("core/geometry.jl")
include("core/mesh.jl")
include("sim/setup.jl")
include("sim/solve.jl")

try
    include("opt/designs.jl")
    include("opt/adjoint.jl")
    include("opt/topopt.jl")
catch e
    @warn "inverse design not available in open source version"
end

include("format.jl")

include("pic/utils.jl")
include("pic/run.jl")

include("plots.jl")
isoval(x::Number) = x
isoval(x::AbstractVector) = mean(x)
isoval(x::AbstractMatrix) = mean(diag(x))
function _run(f, args...; kwargs...)
    try
        f(args...; kwargs...)
        set_job_status("success")
    catch e
        set_job_status("failed")
        println(e)
    end
end

function julia_main()::Cint
    println("running Luminescent FDTD $(VERSION)")
    if isempty(ARGS)
        # println("no arguments provided")
    else
        act, path = ARGS
        if act == "movie"
            _run(makemovie, path)
        else
            _run(picrun, path)
        end
    end
    return 0
end