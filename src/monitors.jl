using LinearAlgebra, UnPack
Base.@kwdef struct Monitor
    span
    k = 0
end
struct MonitorConfig

    idxs
    ki
end
function get(sol::AbstractArray, m,)
    @unpack ki, idxs = m
    [
        begin
            v = map(sol) do v
                v[i][idxs...]
            end
            # isempty(dims) ? vec(v) : reshape(v, dims..., size(v, 2))
        end for i = ki
    ]
end
