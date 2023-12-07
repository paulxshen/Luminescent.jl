using LinearAlgebra, UnPack
struct Monitor
    span
    k
end

# function get(sol::AbstractVector, m, t=eachindex(sol))
#     f = idxs = 0
#     ignore() do
#         @unpack f, idxs = m
#     end
#     r = (;)
#     for f = f
#         r = merge(r, (;
#             f => begin
#                 map(t) do t
#                     sol[t][f][idxs...]
#                 end
#             end
#             # isempty(sz) ? vec(v) : reshape(v, sz..., size(sol, 2))
#         ))
#     end
#     r
# end
function get(sol::AbstractArray, m, t=axes(sol)[end])
    sz = k = i = idxs = 0
    ignore() do
        @unpack sz, k, i, idxs = m
    end
    [
        begin

            v = sol[idxs..., i, t]
            isempty(sz) ? vec(v) : reshape(v, sz..., size(v, 2))
        end for i = i
    ]
end
# function get(sol::AbstractMatrix, m, t)
#     @unpack sz, f, fi, idxs = m
#     r = (;)
#     for f = f
#         r = merge(r, (; f => begin
#             v = sol[fi[f], t]
#         end))
#         # end
#         # r = merge(r, (; f => stack(map(sol.(t)) do v
#         #     v = v[idxs[f]]
#         #     isempty(sz) ? v[1] : reshape(v, sz)
#         # end)))
#     end
#     r
# end
# function get(sol::AbstractVector, m, t)
#     @unpack sz, f, fi, idxs = m
#     r = (;)
#     for f = f
#         r = merge(r, (; f => begin
#             v = sol[fi[f]]
#             isempty(sz) ? vec(v) : reshape(v, sz..., size(sol, 2))
#         end))
#     end
#     r
# end
# function get(sol, m, t)
#     map(m) do m
#         @unpack sz, f, idxs = m
#         r = (;)
#         for f = f
#             r = merge(r, (; f => map(sol.(t)) do v
#                 v = v[f][idxs]
#             end))
#             # r = merge(r, (; f => stack(map(sol.(t)) do v
#             #     v = v[idxs[f]]
#             #     isempty(sz) ? v[1] : reshape(v, sz)
#             # end)))
#         end
#         r
#     end
# end

# function make(m, sol,)
#     map(m) do m
#         @unpack sz, idxs = m
#         r = (;)
#         for f = keys(idxs)
#             r = merge(r, (; f => begin
#                 v = sol[idxs[f], :]
#                 isempty(sz) ? vec(v) : reshape(v, sz..., size(sol, 2))
#             end))
#         end
#         r
#     end
# end