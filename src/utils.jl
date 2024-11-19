getaxisperm(dims) =
    if dims == 1
        [2, 3, 1]
    elseif dims == 2
        [-1, 3, 2]
    elseif dims == 3
        [1, 2, 3]
    end


function permutexyz(d, p)
    _p = @ignore_derivatives invperm(p)
    namedtuple([
        begin
            k = string(k)
            i = findfirst("xyz", k[end])
            Symbol(k[1:end-1] * "xyz"[p[i]]) => sign(p[i]) * permutedims(v, _p)
        end for (k, v) = pairs(d)
    ])
end