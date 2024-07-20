struct PointCloudMonitor <: AbstractMonitor
    p
    n
    w
    c
    s
    sz
    label
end
sphcoords(m::PointCloudMonitor) = m.s
normal(m::AbstractMonitor) = m.n

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
Base.string(m::PointCloudMonitor) =
    """
    $(m.label): point cloud monitor, centered at $(m.c|>d2), $(size(m.p,2)) points"""

struct PointCloudMonitorInstance <: AbstractMonitorInstance
    i
    roi
    n
    w
    c
    inbounds
    label
end
inbounds(m::PointCloudMonitorInstance) = m.inbounds
Base.length(m::PointCloudMonitorInstance) = 1
function MonitorInstance(m::PointCloudMonitor, dx, sz, common_left_pad_amount, flb, fl; F=Float32)
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
    i = p .+ common_left_pad_amount
    roi = NamedTuple([k => p .+ v for (k, v) = pairs(fl)])
    c = round.(c / dx) .+ 1 + common_left_pad_amount
    PointCloudMonitorInstance(i, roi, F.(n), w, c, inbounds, label)
end
function flux(u, m::PointCloudMonitorInstance)
    @unpack n = m
    E = [[field(u, k)[v...] for v = eachcol(m.roi[k])] for k = keys(u[:E])]
    H = [[field(u, k)[v...] for v = eachcol(m.roi[k])] for k = keys(u[:H])]

    dot.(stack(E × H) |> eachrow, eachcol(n))
end
power(u, m::PointCloudMonitorInstance) = flux(u, m) ⋅ m.w
# power(m, u) = sum(sum(pf([u[i...] for (u, i) = zip(u, collect(values(m.i)))]) .* normal(m)))
