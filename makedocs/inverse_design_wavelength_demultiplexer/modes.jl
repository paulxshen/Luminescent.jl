using VectorModesolver
using BSON: @save, @load

ϵ1 = 2.25
ϵ2 = 12.25
ϵbase = ϵclad = 2.25
ϵcore = 12.25
wwg = 0.5
hwg = 0.25
hsub = hclad = lm = 0.25
w = wwg + 2lm
h = hwg + 2lm
dx = 0.05
lb = [-w / 2, -lm]
ub = [w / 2, hwg + lm]
function (ε::εtype)(x, y)

    if (lm <= x <= w - lm) && (lm <= y <= h - lm)
        return (ϵ2, 0.0, 0.0, ϵ2, ϵ2)
    end

    return (ϵ1, 0.0, 0.0, ϵ1, ϵ1)
end

λ1 = 1.55
λ2 = 1
function main(λ)
    ε = εtype()
    x = [i for i in 0:dx:w]
    y = [i for i in 0:dx:h]
    neigs = 1
    tol = 1e-8
    boundary = (0, 0, 0, 0)
    solver = VectorialModesolver(λ, x, y, boundary, ε)
    modes = solve(solver, neigs, tol)

end
modes1_ = main(λ1)
modes2_ = main(λ2)
modes1, modes2 = [[NamedTuple([k => getget(m, k)' for k = [:Ex, :Ey, :Ez, :Hx, :Hy, :Hz]]) for m = m] for m = (modes1_, modes2_)]
@save "$(@__DIR__)/modes.bson" modes1 modes2 λ1 λ2 dx ub lb hsub hwg wwg hclad ϵbase ϵclad ϵcore h w
plot_mode_fields(modes1_[1])
plot_mode_fields(modes2_[1])
