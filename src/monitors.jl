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
    w
    c
    s
    sz
    label
end
sphcoords(m::PointCloudMonitor) = m.s

"""
    function Monitor(c, L; normal=nothing, label="")
    function Monitor(c, lb, ub; normal=nothing, label="")

Constructs monitor which can span a point, line, surface, volume or point cloud monitoring fields or power. 

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

Monitor(p::AbstractMatrix; normals=nothing, weights=nothing, label="") = PointCloudMonitor(p, normals, weights, mean(p, dims=2), nothing, nothing, label)

Monitor(v::AbstractVector; kw...) = Monitor(Matrix(v); kw...)
# Spherical monitor center `c` radius `r` commonly used for antenna pattern . Instantiated to `PointCloudMonitorInstance` consisting of `N` points uniformly sampled on sphere. Points outside simulation domain are automatically discarded 
# function SphereMonitor(c, r; N=256, label="")
"""
    function SphereMonitor(c, r; dθ=5°, dϕ=5°, label="")

Spherical monitor center `c` radius `r` commonly used for antenna pattern . Instantiated to `PointCloudMonitorInstance` . Points outside simulation domain are automatically discarded 
    """
function SphereMonitor(c, r; dθ=2°, dϕ=2°, N=256, nθ=36, nϕ=18, label="")
    # n = stack(sample(Sphere((0, 0, 0), 1), HomogenousSampling(N)))\
    θlims = (0, 360°)
    ϕlims = (0, 180°)
    θ = range(θlims...; step=dθ)
    ϕ = range(ϕlims...; step=dϕ)
    nθ = length(θ)
    nϕ = length(ϕ)
    samples = getproperty.(sample(Sphere((0, 0, 0), 1), RegularSampling(nθ, nϕ)) |> collect, :coords)
    # samples = getproperty.(sample(Sphere((0, 0, 0), 1), HomogeneousSampling(N)) |> collect, :coords)
    # samples=Base.product(θ,ϕ)
    sz = (nθ, nϕ)
    n = stack(samples)
    p = c .+ n * r
    t = SphericalFromCartesian()
    s = stack([[sph.r, sph.θ, sph.ϕ] for sph = t.(samples)])
    # A = 4π * r^2
    # w = A / N
    # w = map(CartesianIndices(sz)) do i, j
    w = map(samples) do (θ, ϕ)
        r^2 * dθ * dϕ * sin(ϕ)
    end
    PointCloudMonitor(p, n, w, c, s, sz, label)
end
Base.string(m::OrthogonalMonitor) =
    """
    $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional monitor, centered at $(m.c|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center, flux normal towards $(m.n|>d2)"""
Base.string(m::PointCloudMonitor) =
    """
    $(m.label): point cloud monitor, centered at $(m.c|>d2), $(size(m.p,2)) points"""

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
Base.size(m::OrthogonalMonitorInstance) = length.(m.i)
support(m::OrthogonalMonitorInstance) = m.v

struct PointCloudMonitorInstance
    i
    fi
    n
    w
    c
    inbounds
    label
end
inbounds(m::PointCloudMonitorInstance) = m.inbounds
Base.length(m::OrthogonalMonitorInstance) = 1
Base.length(m::PointCloudMonitorInstance) = 1

function MonitorInstance(m::OrthogonalMonitor, dx, sz, lc, flb, fl; F=Float32)
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
    i = range.((1 .+ lc + c + lb), (lc + c + ub))
    i = convert(Vector{Any}, i)
    i[singletons] .= first.(getindex.((i,), singletons))
    fi = Dict([k => i .+ v for (k, v) = pairs(flb)])
    # fi = (; [k => i .+ v for (k, v) = pairs(flb)]...)

    n = isnothing(n) ? n : F.(n)
    OrthogonalMonitorInstance(d, i, lc + c, fi, n, F(dx), F(v), m.label)
end

function MonitorInstance(m::PointCloudMonitor, dx, sz, lc, flb, fl; F=Float32)
    @unpack p, n, c, w, label = m
    p = round.(p / dx) .+ 1
    _, i = size(p)
    inbounds = filter(1:i) do i
        v = p[:, i]
        all(v .> 0) && all(v .<= sz)
    end
    p = p[:, inbounds]
    n = n[:, inbounds]
    # p = stack(@. v[1])
    # n = stack(@. v[2])
    i = p .+ lc
    fi = NamedTuple([k => p .+ v for (k, v) = pairs(fl)])
    c = round.(c / dx) .+ 1 + lc
    PointCloudMonitorInstance(i, fi, F.(n), w, c, inbounds, label)
end

# function Base.getindex(a, m::OrthogonalMonitorInstance, )
#     a[m.fi[k]...]
# end

function field(u, k)
    # i, j = fij[k]
    # u[i][j]
    if k in keys(u)
        u[k]
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
    else
        u[k[1]][k]
    end
    # error("invalid field")
end
# if k in keys(u)
#     u[k]

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
    r = _get(u, k, m)
    if isnothing(r)
        return field(u, k)[m.fi[k]...]
    end
    r
    # field(u, k, m)
end
function field(u, k, m::PointCloudMonitorInstance,)
    r = _get(u, k, m)
    if isnothing(r)
        r = [field(u, k)[v...] for v = eachcol(m.fi[k])]
    end
    # if !isnothing(m.sz)
    #     r = reshape(r, m.sz)
    # end
    r
end
function _get(u, k, m,)
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

# get
# Base.getindex(a::Union{AbstractArray,GPUArraysCore.AbstractGPUArray}, m::MonitorInstance) = a[m.i...]
# Base.getindex(a::GPUArraysCore.AbstractGPUArray, m::MonitorInstance) = a[m.i...]
# Base.getindex(a, m::MonitorInstance, k) = a[m.fi[k]...]

# power(m, u) = sum(sum(pf([u[m.i[k]...] for (u, k) = zip(u, fk)]) .* m.n))
"""
    function flux(m::MonitorInstance, u)

 Poynting flux profile passing thru monitor 
"""
function flux(u, m::OrthogonalMonitorInstance)
    d = ndims(m)
    E = [u[:E][k][m.fi[k]...] for k = keys(u[:E])]
    H = [u[:H][k][m.fi[k]...] for k = keys(u[:H])]

    sum((E × H) .* m.n)
end
function flux(u, m::PointCloudMonitorInstance)
    @unpack n = m
    E = [[field(u, k)[v...] for v = eachcol(m.fi[k])] for k = keys(u[:E])]
    H = [[field(u, k)[v...] for v = eachcol(m.fi[k])] for k = keys(u[:H])]

    dot.(stack(E × H) |> eachrow, eachcol(n))
end

"""
    function power(m::MonitorInstance, u)

total power (Poynting flux) passing thru monitor surface
"""
function power(u, m::OrthogonalMonitorInstance)
    # @assert ndims(m) == 2
    m.v * mean(flux(u, m))
end
power(u, m::PointCloudMonitorInstance) = flux(u, m) ⋅ m.w
# power(m, u) = sum(sum(pf([u[i...] for (u, i) = zip(u, collect(values(m.i)))]) .* m.n))

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