module Luminescent
include("main.jl")
export sandwich
export Periodic, PML, PEC, PMC, InPad, OutPad
export PlaneWave, GaussianBeam, Source
export Monitor, SphereMonitor, power, flux, field, support, sphcoords, inbounds
export setup, apply, apply!, °, update!, update
export calibrate_mode, collapse_mode
export port_number, mode_number, shorten_key, simplify_sparams, sparam_family
export solve, dispersion_compensation
end

# ]add UnPack,StatsBase, DataStructures, LinearAlgebra, Random,  Interpolations, Flux,Zygote, Optim,ArrayPadding, Jello