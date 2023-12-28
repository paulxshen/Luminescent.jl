module fdtd_prerelease
using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Flux, Optim, DataStructures
using Flux: withgradient
using Optim: Options, minimizer
using Zygote: ignore, Buffer
using Jello, ArrayPadding

F = Float32
Random.seed!(1)
include("utils.jl")
include("del.jl")
include("maxwell.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("fdtd.jl")

export reindex
export Del, stepTMz
export Periodic, PML, PEC, PMC, Padding
export PlaneWave, GaussianBeam, UniformSource, CenteredSource, place
export Monitor
export setup, apply
end

# ]add UnPack,StatsBase, DataStructures, LinearAlgebra, Random,  Interpolations, Flux,Zygote, Optim,ArrayPadding, Jello