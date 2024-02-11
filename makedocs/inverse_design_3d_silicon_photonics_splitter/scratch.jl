error()
# @show x|>re|>metrics

# train surrogate
Random.seed!(1)
runs = []
i = 1
function f!(x)
    @time begin
        model = re(x)
        mp = metrics(model)
        l = loss(mp)
        global i
        # if i == 1
        push!(runs, [x, mp, l])
        # else
        #     runs[i] .= [x0, mp, l]
        # end
        # i +=1
        l
    end
end

iterations = 24
f!(x)
# runs = similar(runs, iterations)

@showtime res = optimize(f!, x0, ParticleSwarm(;
        n_particles=16), Optim.Options(; f_tol=0, iterations, show_every=1, show_trace=true))
xgf = minimizer(res)
x = deepcopy(xgf)

error()
X = Base.stack(getindex.(runs, 1))
Y = Base.stack(getindex.(runs, 2))

nl = leakyrelu
n = size(X, 1)
m, N = size(Y)
nn = Chain(Dense(n, 2n, nl), Dense(2n, 4n, nl), Dense(4n, m))

opt = Adam(0.1)
opt_state = Flux.setup(opt, nn)
n = 100
for i = 1:n
    l, (dldm,) = withgradient(nn -> Flux.mae(Y, nn(X)), nn)
    Flux.update!(opt_state, nn, dldm)
    (i % 25 == 0 || i == n) && println("$i $l")
end

opt = Adam(0.1)
opt_state = Flux.setup(opt, x)
n = 100
@show f(x)
for i = 1:n
    l, (dldm,) = withgradient(x -> mae(
            nn(x)[2:end],
            [tp / 2, tp / 2, -tp, tp, zeros(F, 4)...]
        ), x)
    Flux.update!(opt_state, x, dldm)
    (i % 25 == 0 || i == n) && println("$i $l")
end
@show loss(nn(x))
@show f(x)
xsg = x
x = deepcopy(xsg)
