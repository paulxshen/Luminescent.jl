# Base.convert(::Type{Float64}, x::ComplexF64) = real(x)
using VectorModesolver

function ε(x::Float64, y::Float64)

    if (2 < x < 5) && (2 < y < 3)
        return (4.0, 0.0, 0.0, 4.0, 4.0)
    end

    return (1.0, 0.0, 0.0, 1.0, 1.0)
end

function main()
    λ = 1.55
    x = [i for i in 0:0.05:7]
    y = [i for i in 0:0.05:5]
    neigs = 3
    tol = 1e-8
    boundary = (0, 0, 0, 0)
    solver = VectorialModesolver(λ, x, y, boundary, ε)
    modes = VectorModesolver.solve(solver, neigs, tol)

    plot_mode_fields(modes[1]) |> display
    plot_mode_fields(modes[2]) |> display
end

main()
