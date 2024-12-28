function cdiff(a; dims=1)
    n = size(a, dims)
    selectdim(a, dims, 3:n) - selectdim(a, dims, 1:n-2)
end

function diffpad(a, vl, vr=vl; dims=1, diff=diff)
    # @assert all(!isnan, a)

    sel = 1:ndims(a) .== dims
    l = !isnothing(vl)
    r = !isnothing(vr)

    sz = Tuple(size(a) + (l + r - 1) * sel)
    b = Buffer(a, sz)
    b[range.(l * sel + 1, sz - r * sel)...] = diff(a; dims)
    pad!(b, vl, l * sel, 0)
    pad!(b, vr, 0, r * sel)
    copy(b)
end

struct Del
    diff
    deltas
    padvals
    function Del(deltas=1, padvals=nothing; diff=diff)
        new(diff, deltas, padvals)
    end
end

function LinearAlgebra.cross(m::Del, d::Map)
    @unpack diff, deltas, padvals = m
    ks = keys(d)
    Δs = @ignore_derivatives [[d(k) for k in ks] for d = deltas]
    ps = @ignore_derivatives [[d(k) for k in ks] for d = padvals]
    as = values(d)
    delcross(diff, Δs, ps, as)
end

function LinearAlgebra.cross(m::Del, v)
    @unpack diff, deltas, padvals = m
    Δs = values.(deltas)
    ps = values.(padvals)
    delcross(diff, Δs, ps, v)
end

function delcross(diff, Δs, ps, as)
    _diffpad(a, p, dims) = diffpad(a, p...; dims, diff)
    N = ndims(as(1))
    if N == 2
        dx, dy = [Δs(i) for i = 1:2]
        if length(as) == 1
            u, = as
            return [_diffpad(u, ps(2)(1), 2) ./ dy(1), -_diffpad(u, ps(1)(1), 1) ./ dx(1)]
        elseif length(as) == 2
            u, v = as
            return [_diffpad(v, ps(1)(2), 1) ./ dx(2) - _diffpad(u, ps(2)(1), 2) ./ dy(1)]
        else
            error()
        end
    elseif N == 3
        dx, dy, dz = [Δs(i) for i = 1:3]
        u, v, w = as
        uy = _diffpad(u, ps(2)(1), 2) ./ dy(1)
        uz = _diffpad(u, ps(3)(1), 3) ./ dz(1)
        vx = _diffpad(v, ps(1)(2), 1) ./ dx(2)
        vz = _diffpad(v, ps(3)(2), 3) ./ dz(2)
        wx = _diffpad(w, ps(1)(3), 1) ./ dx(3)
        wy = _diffpad(w, ps(2)(3), 2) ./ dy(3)
        return [wy - vz, uz - wx, vx - uy]
    end
    error()
end