using LinearAlgebra, UnPack
struct Monitor
    span
    fn
end

function make(m, sol,)
    map(m) do m
        @unpack sz, idxs = m
        r = (;)
        for f = keys(idxs)
            r = merge(r, (; f => begin
                v = sol[idxs[f], :]
                isempty(sz) ? vec(v) : reshape(v, sz..., size(sol, 2))
            end))
        end
        r
    end
end