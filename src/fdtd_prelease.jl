module fdtd_prelease
using UnPack, LinearAlgebra, Random, Images, Interpolations, Flux, Optim, DataStructures
using Flux: withgradient
using Optim: Options, minimizer
using Zygote: ignore, Buffer
using BSON: @save, @load
using Jello, ArrayPadding

F = Float32
Random.seed!(1)
include("startup.jl")
include("del.jl")
include("maxwell.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("fdtd.jl")

export Del, stepTMz
export Periodic, PML, PEC, PMC
export PlaneWave, GaussianBeam, UniformSource
export Monitor
export setup
end
