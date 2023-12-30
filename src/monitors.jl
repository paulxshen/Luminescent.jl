using LinearAlgebra, UnPack
"""
    function Monitor(span, normal=nothing)

Constructs monitor which can span a point, line, surface, or volume

Args
- span
- normal: flux monitor direction
"""
struct Monitor
    span
    normal
    function Monitor(span, normal=nothing)
        new(span, normal)
    end
end
# struct MonitorConfig

#     idxs
#     ki
# end
# function get(sol::AbstractArray, m,)
#     @unpack ki, idxs = m
#     [
#         begin
#             v = map(sol) do v
#                 v[i][idxs...]
#             end
#             # isempty(dims) ? vec(v) : reshape(v, dims..., size(v, 2))
#         end for i = ki
#     ]
# end
