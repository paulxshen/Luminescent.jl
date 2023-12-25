using NPZ
using BSON: @save, @load
F = Float32
function invrs_load(fn)
    f(x) = 0.1 < x < 0.9
    base = F.(npzread("examples/" * fn * ".npy"))
    _c = 0
    base = base[1+_c:end-_c, 1+_c:end-_c]
    design_start = Tuple(findfirst(f, base))
    design_sz = 1 .+ Tuple(Tuple(findlast(f, base)) .- design_start)
    dx = 1.6f0 / 40
    base[[a:b for (a, b) = zip(design_start, design_start .+ design_sz .- 1)]...] .= 0

    # 
    @save "examples/$(fn).bson" base design_start design_sz dx
end
for fn = ["waveguide_bend"]
    invrs_load(fn)
end