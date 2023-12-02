using NPZ
f(x) = 0.1 < x < 0.9
F = Float32
function invrs_load(fn)
    base = F.(npzread(fn * ".npy"))
    mstart = Tuple(findfirst(f, base))
    msz = 1 .+ Tuple(Tuple(findlast(f, base)) .- mstart)
    dx_ = 1.6f0 / 40

    base, mstart, msz, dx_
end