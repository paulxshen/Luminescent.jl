abstract type AbstractMonitor end
struct Monitor <: AbstractMonitor
    c
    lb
    ub
    n
    tangent
    label

    wavelength_modes
    meta
end

wavelengths(m::Monitor) = keys(m.wavelength_modes)

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
function Monitor(c, L, normal=nothing; label="", kw...)
    Monitor(c, -L / 2, L / 2, normal, label, nothing, kw)
end

function ModalMonitor(wavelength_modes::Dictlike, c, normal, tangent, L; wavelength=1, wavelengths=[wavelength], label="", kw...)
    Monitor(c, -L / 2, L / 2, normal, tangent, label, wavelength_modes, kw)
end

Base.string(m::Monitor) =
    """
    $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional monitor, centered at $(m.c|>d2), physical size of $((m.ub-m.lb)|>d2) relative to center, flux normal towards $(normal(m)|>d2)"""

abstract type AbstractMonitorInstance end
mutable struct MonitorInstance <: AbstractMonitorInstance
    d
    roi
    frame
    dx
    v
    center
    label

    wavelength_modes
end
@functor MonitorInstance (wavelength_modes,)
Base.ndims(m::MonitorInstance) = m.d
# Base.size(m::MonitorInstance) = length.(m.i) |> Tuple
area(m::MonitorInstance) = m.v
wavelengths(m::MonitorInstance) = Porcupine.keys(m.wavelength_modes)


Base.length(m::MonitorInstance) = 1

frame(m::MonitorInstance) = m.frame
normal(m::MonitorInstance) = frame(m)[3][1:length(m.center)]

function MonitorInstance(m::Monitor, dx, field_origin, common_left_pad_amount, sz, ; F=Float32)
    @unpack n, lb, ub, c, tangent, wavelength_modes, = m
    n, lb, ub, c, tangent = F.((n, lb, ub, c, tangent))
    L = ub - lb
    D = length(c)
    d = length(lb)
    A = isempty(L) ? 0 : prod(L,)
    if D == 2
        frame = [[-n[2], n[1], 0], [0, 0, 1], [n..., 0]]

    elseif D == 3
        frame = [tangent, n × tangent, n]
    end
    lb = sum(lb .* frame[1:d])[1:D]
    ub = sum(ub .* frame[1:d])[1:D]
    v = zip(lb, ub)
    lb = minimum.(v)
    ub = maximum.(v)
    L = ub - lb

    roi = dict([k => begin
        a = (c + lb - o) / dx + 1.5

        [l == 0 ? a : (a + (0:sign(l):l-sign(l))) for (a, l) in zip(a, round(L / dx))]
        # p = (c + lb - o) / dx + 1.5
        # i = floor(p)
        # w = 1 - mean(abs.(p - i))
        # w = (w, 1 - w)
        # i = [i .+ map(round(L / dx)) do n
        #     n == 0 ? 0 : (0:sign(n):n-sign(n))
        # end for i = (i, i + 1)]
        # [(F(w), i) for (w, i) in zip(w, i) if w > 0]
    end for (k, o) = pairs(field_origin)])
    n = isnothing(n) ? n : F.(n |> normalize)
    _center = round(c / dx) + 1 + common_left_pad_amount

    MonitorInstance(d, roi, frame, F(dx), F(A), _center, m.label, F(wavelength_modes))
end


"""
    function field(u, k, m)
    
queries field, optionally at monitor instance `m`

Args
- `u`: state
- `k`: symbol or str of Ex, Ey, Ez, Hx, Hy, Hz, |E|, |E|2, |H|, |H|2
- `m`
"""
function field(a::AbstractArray, k, m)
    # sum([w * a[i...] for (w, i) = m.roi[k]])
    getindexf(a, m.roi[k]...)
end

function field(u::Dictlike, k, m)
    field(u(k), k, m)
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
#         #     return sum([w * a[i...] for (w, i) = m.roi[k]])
#         # elseif isa(m, PointCloudMonitorInstance)
#         #     return [a[v...] for v = eachcol(m.roi[k])]
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