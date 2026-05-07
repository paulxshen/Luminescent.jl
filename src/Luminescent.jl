module Luminescent
include("main.jl")
export _run, picrun, makemovie, set_job_status, set_server_status
# export Periodic, PML, PEC, PMC, InPad, OutPad
# export PlaneWave, GaussianBeam, Source, Source, keepxy
# export AbstractMonitor, AbstractMonitor, SphereMonitor, power, flux, field, support, sphcoords, inbounds
# export setup, apply, apply!, °, update!, update
# export calibrate_mode, collapse_mode
# export port_number, mode_number, shorten_key, simplify_sparams, sparam_family
# export solve, dispersion_compensation, quickie
# export gfrun, lastrun, cpu, gpu
end

# ]add UnPack,StatsBase, DataStructures, LinearAlgebra, Random,  Interpolations, Flux,Zygote, Optim,ArrayPadding, Jello