using VectorModesolver

function (ε::εtype)(x::Float64, y::Float64)

    if (0.75 < x < 1.75) && (0.95 < y < 1.55)
        return (4.0, 0.0, 0.0, 4.0, 4.0)
    end

    return (1.0, 0.0, 0.0, 1.0, 1.0)
end

function main()
    ε = εtype()
    λ = 1.55
    x = [i for i in 0:0.03:2.5]
    y = [i for i in 0:0.05:2.5]
    neigs = 1
    tol = 1e-8
    boundary = (0, 0, 0, 0)
    solver = VectorialModesolver(λ, x, y, boundary, ε)
    modes = solve(solver, neigs, tol)

end

modes = main()
plot_mode_fields(modes[1])
