module Luminescent
include("main.jl")
export sandwich
export maxwell_update!, maxwell_update!, maxwell_update
# stepTM, step1, , stepTE
export Periodic, PML, PEC, PMC, InPad, OutPad
export PlaneWave, GaussianBeam, Source, place, place!
export Monitor, power, flux, field, support, sphcoords, inbounds
export maxwell_setup, apply, apply!, °
end

# ]add UnPack,StatsBase, DataStructures, LinearAlgebra, Random,  Interpolations, Flux,Zygote, Optim,ArrayPadding, Jello