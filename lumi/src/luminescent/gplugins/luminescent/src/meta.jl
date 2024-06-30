using DiffEqFlux, DifferentialEquations, ArrayPadding
include("del.jl")

dx = 1 / 16.0f0
lx = ly = lz = 4.0f0
dims = round(Int, [lx, ly, lz] / dx)

ϵ = ones(dims)
μ = ones(dims)
σ = zeros(dims)

u0 = zeros(dims..., 6)
τ = 1.0f0
f(t) = exp(-((t - τ))^2)
datasize = 8
tspan = (0.0f0, 6.0f0)
tsteps = range(tspan[1], tspan[2], length=datasize)

function trueODEfunc(du, u, p, t)
    Ex, Ey, Ez, Hx, Hy, Hz = eachslice(u)

    Ex = cat(zeros(dims[2:3]), Ex, dims=3)
    Ey = cat(zeros(dims[2:3]), Ey, dims=3)
    Ez = cat(f(t) * ones(dims[2:3]), Ez, dims=3)

    Hx = cat(zeros(dims[2:3]), Hx, dims=3)
    Hy = cat(zeros(dims[2:3]), Hy, dims=3)
    Hz = cat(zeros(dims[2:3]), Hz, dims=3)

    Ex = pad(Ex, :periodic, (0, 1, 1), (0, 0, 0))
    Ey = pad(Ey, :periodic, (0, 1, 1), (0, 0, 0))
    Ez = pad(Ez, :periodic, (0, 1, 1), (0, 0, 0))
    Hx = pad(Hx, :periodic, (0, 0, 0), (0, 1, 1),)
    Hy = pad(Hy, :periodic, (0, 0, 0), (0, 1, 1),)
    Hz = pad(Hz, :periodic, (0, 0, 0), (0, 1, 1),)


    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u .^ 3)'true_A)'
end

prob_trueode = ODEProblem(trueODEfunc, u0, tspan)
ode_data = Array(solve(prob_trueode, Tsit5(), saveat=tsteps))