function vec2cmat(v, dims)
    a = reshape(v, dims..., 2)
    complex.(a[:, :, 1], a[:, :, 2])
end
function cmat2vec(m)
    v = vec(m)
    vcat(real(v), imag(v))
end
function f(x)
    model.a .= vec2cmat(x, size(model.a))
end

function g!(storage, x)
    storage .= cmat2vec(g)
end
function fg!(storage, x)
    storage .= cmat2vec(g)
    l
end
x0 = cmat2vec(model.a)
od = OnceDifferentiable(f, g!, fg!, x0)
minimizer(optimize(od, x0, LBFGS(), iterations=nepochs))