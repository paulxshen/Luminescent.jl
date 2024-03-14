# For example, `a::Array[m::MonitorInstance, :Ex]` queries Ex on monitor m.
fij = (
    Ex=(1, 1),
    Ey=(1, 2),
    Ez=(1, 3),
    Hx=(2, 1),
    Hy=(2, 2),
    Hz=(2, 3),
)
struct OrthogonalMonitor
    c
    lb
    ub
    n
    label
end

struct PointCloudMonitor
    p
    n
    c
    label
end

"""
    function Monitor(c, L; normal=nothing, label="")
    function Monitor(c, lb, ub; normal=nothing, label="")

Constructs monitor which can span a point, line, surface, volume or point cloud monitoring fields or power_flux. 

Args
- c: origin or center of monitor
- L: physical dimensions of monitor
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c

- normal: flux monitor direction (eg normal to flux surface)
"""
function Monitor(c, L; normal=nothing, label="")
    OrthogonalMonitor(c, -L / 2, L / 2, normal, label)
end

function Monitor(c, lb, ub; normal=nothing, label="")
    OrthogonalMonitor(c, lb, ub, normal, label)
end

Monitor(p::AbstractMatrix; normal=nothing, label="") = PointCloudMonitor(p, normal, mean(p, dims=2), label)

Monitor(v::AbstractVector; kw...) = Monitor(Matrix(v); kw...)
Base.string(m::OrthogonalMonitor) =
    """
    $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional monitor, centered at $(m.c|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center, flux normal towards $(m.n|>d2)"""

struct OrthogonalMonitorInstance
    d
    i
    c
    fi
    n
    dx
    v
    label
end
# @functor OrthogonalMonitorInstance (n,)
Base.ndims(m::OrthogonalMonitorInstance) = m.d

struct PointCloudMonitorInstance
    i
    fi
    n
    c
    label
end

function MonitorInstance(m::OrthogonalMonitor, dx, lc, flb, fl; F=Float32)
    @unpack n, = m
    L = m.ub - m.lb
    singletons = findall(L .≈ 0,)
    d = length(L) - length(singletons)
    deleteat!(L, singletons)
    v = isempty(L) ? 0 : prod(L,)
    c = round.(Int, m.c / dx)
    lb = round.(Int, m.lb / dx)
    ub = round.(Int, m.ub / dx)

    # idxs = Dict([
    #     k => map(sizes[k],  c, lb, ub) do s, c, lb, ub
    #         max(1, lc + c + lb):min(s, lc + c + ub)
    #     end for k = fk
    # ])
    i = range.((1 .+ lc + c + lb), (1 .+ lc + c + ub))
    i = convert(Vector{Any}, i)
    i[singletons] .= first.(getindex.((i,), singletons))
    fi = Dict([k => i .+ v for (k, v) = pairs(flb)])
    # fi = (; [k => i .+ v for (k, v) = pairs(flb)]...)

    n = isnothing(n) ? n : F.(n)
    OrthogonalMonitorInstance(d, i, lc + c, fi, n, F(dx), F(v), m.label)
end

function MonitorInstance(m::PointCloudMonitor, dx, lc, flb, fl; F=Float32)
    @unpack p, n, c, label = m
    p = round.(p / dx) .+ 1
    i = p .+ lc
    fi = NamedTuple([k => p .+ v for (k, v) = pairs(fl)])
    c = round.(c / dx) .+ 1 + lc
    PointCloudMonitorInstance(i, fi, F.(n), c, label)
end

# function Base.getindex(a, m::OrthogonalMonitorInstance, )
#     a[m.fi[k]...]
# end

function field(u, k)
    if k in keys(fij)
        u[k[1]][k]
        # i, j = fij[k]
        # u[i][j]
    elseif k == "|E|2"
        sum(u[1]) do a
            a .^ 2
        end
    elseif k == "|E|"
        sqrt.(field(u, "|E|2"))
    elseif k == "|H|2"
        sum(u[2]) do a
            a .^ 2
        end
    elseif k == "|H|"
        sqrt.(field(u, "|H|2"))
    end
end

"""
    function field(u, k)
    function field(u, k, m)

queries field, optionally at monitor instance `m`

Args
- `u`: state
- `k`: symbol or str of Ex, Ey, Ez, Hx, Hy, Hz, |E|, |E|2, |H|, |H|2
- `m`
"""
function field(u, k, m::OrthogonalMonitorInstance,)
    if k in keys(fij)
        i, j = fij[k]
        return u[i][j][m.fi[k]...]
    end
    field(u, k, m)
end
function field(u, k, m::PointCloudMonitorInstance,)
    if k in keys(fij)
        i, j = fij[k]
        return [u[i][j][v...] for v = eachcol(m.fi[k])]
    end
    field(u, k, m)
end
function field(u, k, m,)
    if k == "|E|2"
        sum(field.(u, (:Ex, :Ey, :Ez), (m,))) do a
            a .^ 2
        end
    elseif k == "|E|"
        sqrt.(field(u, "|E|2", m))
    elseif k == "|H|2"
        sum(field.(u, (:Hx, :Hy, :Hz), (m,))) do a
            a .^ 2
        end
    elseif k == "|H|"
        sqrt.(field(u, "|H|2", m))
    end
end
# Base.getindex(a::Union{AbstractArray,GPUArraysCore.AbstractGPUArray}, m::MonitorInstance) = a[m.i...]
# Base.getindex(a::GPUArraysCore.AbstractGPUArray, m::MonitorInstance) = a[m.i...]
# Base.getindex(a, m::MonitorInstance, k) = a[m.fi[k]...]

# power_flux(m, u) = sum(sum(pf([u[m.i[k]...] for (u, k) = zip(u, fk)]) .* m.n))
"""
    function power_flux_density(m::MonitorInstance, u)

 power_flux density (avg Poynting flux) passing thru monitor 
"""
function power_flux_density(m::OrthogonalMonitorInstance, u)
    d = ndims(m)
    E = [u[:E][k][m.fi[k]...] for k = keys(u[:E])]
    H = [u[:H][k][m.fi[k]...] for k = keys(u[:H])]

    r = mean(sum((E × H) .* m.n))
    r
end

"""
    function power_flux(m::MonitorInstance, u)

total power_flux (Poynting flux) passing thru monitor surface
"""
function power_flux(m::OrthogonalMonitorInstance, u)
    @assert ndims(m) == 2
    m.v * power_flux_density(m, u)
end
# power_flux(m, u) = sum(sum(pf([u[i...] for (u, i) = zip(u, collect(values(m.i)))]) .* m.n))

function monitors_on_box(c, L)
    ox, oy, oz = c
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