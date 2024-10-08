using VectorModesolver
using BSON: @save, @load

ϵ1 = 2.25
ϵ2 = 12.25
ϵbase = ϵclad = 2.25
ϵcore = 12.25
wwg = 0.5
hwg = 0.25
hbase = hclad = lm = hm = 0.25
wm = 0.5
w = wwg + 2wm
h = hwg + 2lm
dx = 0.05
lb = [-w / 2, -hm]
ub = [w / 2, hwg + hm]
function (ε::εtype)(x, y)

    if (lm <= x <= w - lm) && (lm <= y <= h - lm)
        return (ϵ2, 0.0, 0.0, ϵ2, ϵ2)
    end

    return (ϵ1, 0.0, 0.0, ϵ1, ϵ1)
end

λ = 1.55
function main()
    ε = εtype()
    x = [i for i in 0:dx:w]
    y = [i for i in 0:dx:h]
    neigs = 1
    tol = 1e-8
    boundary = (0, 0, 0, 0)
    solver = VectorialModesolver(λ, x, y, boundary, ε)
    modes = solve(solver, neigs, tol)

end
modes_ = main()
modes = [NamedTuple([k => transpose(getfield(m, k),) for k = [:neff, :Ex, :Ey, :Ez, :Hx, :Hy, :Hz]]) for m = modes_]
@save "$(@__DIR__)/modes.bson" modes λ dx ub lb hbase hwg wwg hclad ϵbase ϵclad ϵcore h w
plot_mode_fields(modes_[1])
