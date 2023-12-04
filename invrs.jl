using NPZ
using BSON: @save, @load
F = Float32
function invrs_load(fn)
    f(x) = 0.1 < x < 0.9
    base = F.(npzread(fn * ".npy"))
    design_start = Tuple(findfirst(f, base))
    design_sz = 1 .+ Tuple(Tuple(findlast(f, base)) .- design_start)
    dx = 1.6f0 / 40
    base[[a:b for (a, b) = zip(design_start, design_start .+ design_sz .- 1)]...] .= 0


    @save "$fn.bson" base design_start design_sz dx
end
for f = ["bend"]
    invrs_load(f)
end