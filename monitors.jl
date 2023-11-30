using LinearAlgebra, UnPack
struct Monitor
    span
    fn
end

function make(m, sol, t)
    map(m) do m
        @unpack sz, idxs = m
        # NamedTuple([f => stack(
        #     map(sol.(t)) do v
        #         v = v[idxs[f]]
        #         if isempty(sz)
        #             v[1]
        #         else
        #             reshape(v, sz)
        #         end
        #     end
        # )
        #             for f = keys(idxs)])

        #         end
        r = (;)
        for f = keys(idxs)
            r = merge(r, (; f => map(sol.(t)) do v
                v = v[idxs[f]][1]
            end))
            # r = merge(r, (; f => stack(map(sol.(t)) do v
            #     v = v[idxs[f]]
            #     isempty(sz) ? v[1] : reshape(v, sz)
            # end)))
        end
        r
    end
end