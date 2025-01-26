abstract type AbstractMonitor end

struct Monitor <: AbstractMonitor
    λmodenums
    λsmode
    λmodes

    center
    L
    center3
    L3
    mask


    dimsperm
    frame

    approx_2D_mode
    # tangent

    tags
    # function Monitor(center, L, normal=nothing; tags...)
    #     new(center, -L / 2, L / 2, normal, nothing, tags)
    # end
    # function Monitor(λmodes::Map, center, normal, tangent, L; tags...)
    #     Monitor(center, -L / 2, L / 2, normal, tangent, λmodes, tags)
    # end
    # function Monitor(a...)
    #     new(a...)
    # end
end

Monitor(center, L, dimsperm, center3=center, L3=L, approx_2D_mode=nothing; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...) =
    Monitor(λmodenums, λsmode, λmodes, center, L, center3, L3, nothing, dimsperm, nothing, approx_2D_mode, tags)

Monitor(mask, frame; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...) =
    Monitor(λmodenums, λsmode, λmodes, nothing, nothing, nothing, nothing, mask, nothing, frame, nothing, tags)

function Base.ndims(m::Monitor)
    v = m.center
    if !isnothing(v)
        length(v)
    else
        ndims(m.mask)
    end
end

struct PlaneMonitor <: AbstractMonitor
    dims
    q

    λmodes
    # meta
    tags
    function PlaneMonitor(q, λmodes; dims, tags...)
        new(dims, q, λmodes, tags)
    end
end



wavelengths(m::Monitor) = keys(m.λmodes)

abstract type AbstractMonitorInstance end
struct MonitorInstance <: AbstractMonitorInstance
    inds
    masks

    frame
    dimsperm
    deltas
    center
    λmodes
    _λmodes
    tags
end
@functor MonitorInstance (λmodes, _λmodes, masks, deltas)
Base.ndims(m::MonitorInstance) = length(m.center)
area(m::MonitorInstance) = m.v
wavelengths(m::MonitorInstance) = keys(m.λmodes)
Base.length(m::MonitorInstance) = 1
frame(m::MonitorInstance) = m.frame
normal(m::MonitorInstance) = frame(m)[3][1:length(m.center)]

function MonitorInstance(m::Monitor, g, ϵ, TEMP, mode_solutions=nothing)
    λmodes, _λmodes, inds, masks, labelpos, = _get_λmodes(m, ϵ, TEMP, mode_solutions, g)

    MonitorInstance(inds, masks, m.frame, m.dimsperm, g.deltas, labelpos, λmodes, _λmodes, m.tags)
end

function MonitorInstance(m::PlaneMonitor, g)
    @unpack L = g
    @unpack dims, q, λmodes, tags = m
    MonitorInstance(Monitor(λmodes, L / 2, -L / 2, L / 2, getdimsperm(dims), tags), g)
end


# """
#     function field(u, k, m)

# queries field, optionally at monitor instance `m`

# Args
# - `u`: state
# - `k`: symbol or str of Ex, Ey, Ez, Hx, Hy, Hz, |E|, |E|2, |H|, |H|2
# - `m`
# """
# function field(a::AbstractArray, k, m)

#     # permutedims(getindexf(a, m.inds[k]...), p)
# end

function field(u::Map, k, m)
    @unpack inds, masks = m
    @nograd k, m, inds, masks
    a = u(k)
    if !isnothing(masks)
        a .* masks(k)
    else
        getindexf(a, m.inds[k]...)#; approx=true)
    end
end

#     if k == "|E|2"
#         sum(field.(u, (:Ex, :Ey, :Ez), (m,))) do a
#             a .|> abs2
#         end
#     elseif k == "|E|"
#         sqrt.(field(u, "|E|2", m))
#     elseif k == "|H|2"
#         sum(field.(u, (:Hx, :Hy, :Hz), (m,))) do a
#             a .|> abs2
#         end
#     elseif k == "|H|"
#         sqrt.(field(u, "|H|2", m))
#     else
#         # a =field( u(k),k,m)
#         # a = u[k[1]][k]
#         # if isnothing(m)
#         #     return a
#         # elseif isa(m, MonitorInstance)
#         #     return sum([w * a[i...] for (w, i) = m.inds[k]])
#         # elseif isa(m, PointCloudMonitorInstance)
#         #     return [a[v...] for v = eachcol(m.inds[k])]
#         # end
#     end
# end

Base.getindex(a::AbstractDict, k, m::MonitorInstance) = field(a, k, m)
Base.getindex(a::NamedTuple, k, m::MonitorInstance) = field(a, k, m)

"""
    function flux(m::MonitorInstance, u)

 Poynting flux profile passing thru monitor 
"""
function flux(u, polarization=nothing)
    E = group(u, :E)
    H = group(u, :H)

    P_TE = E.Ex .* conj.(H.Hy)
    P_TM = -E.Ey .* conj.(H.Hx)

    if polarization == :TE
        return P_TE
    elseif polarization == :TM
        return P_TM
    else
        return P_TE + P_TM
    end
    #     # (E × H) ⋅ normal(m)
    #    return sum((E × conj.(H)) .* normal(m))
end

"""
    function power(m::MonitorInstance, u)

total power (Poynting flux) passing thru monitor surface
"""
function power(u, m::MonitorInstance)
    # @assert ndims(m) == 2
    m.v * mean(flux(u, m))
end

function monitors_on_box(center, L)
    ox, oy, oz = center
    lx, ly, lz = L
    rx, ry, rz = L / 2
    [
        Monitor([ox - rx, oy, oz], [0, ly, lz], [-1, 0, 0]),
        Monitor([ox + rx, oy, oz], [0, ly, lz], [1, 0, 0]),
        Monitor([ox, oy - ry, oz], [lx, 0, lz], [0, -1, 0]),
        Monitor([ox, oy + ry, oz], [lx, 0, lz], [0, 1, 0]),
        Monitor([ox, oy, oz - rz], [lx, ly, 0,], [0, 0, -1,]),
        Monitor([ox, oy, oz + rz], [lx, ly, 0,], [0, 0, 1,]),
    ]
end