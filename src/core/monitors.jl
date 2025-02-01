abstract type AbstractMonitor end

mutable struct Monitor <: AbstractMonitor
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
mutable struct MonitorInstance <: AbstractMonitorInstance
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
Base.getindex(a::AbstractDict, k, m::MonitorInstance) = field(a, k, m)
Base.getindex(a::NamedTuple, k, m::MonitorInstance) = field(a, k, m)

