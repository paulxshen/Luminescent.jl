abstract type AbstractMonitor end
@getr AbstractMonitor

struct BaseMonitor <: AbstractMonitor
    name
    origin
    meta
end
Base.ndims(m::AbstractMonitor) = length(m.origin) - 1

mutable struct PlaneMonitor <: AbstractMonitor
    base_monitor::BaseMonitor
    start
    stop
    frame
    voltage_line
    current_loop
    function PlaneMonitor(name, origin, start, stop, frame; voltage_line=nothing, current_loop=nothing, meta...)
        new(BaseMonitor(name, origin, meta), start, stop, frame, voltage_line, current_loop)
    end
end

mutable struct PointsMonitor <: AbstractMonitor
    base_monitor::BaseMonitor
    function PointsMonitor(name, origin; meta...)
        new(BaseMonitor(name, origin, meta))
    end
end

mutable struct SphereMonitor <: AbstractMonitor
    points_monitor::PointsMonitor
    frame
    radius
    θ
    φ
    function SphereMonitor(name, origin, radius::Real, frame=Matrix(identity, 3, 3); angres=deg2rad(10), θmin=1f-3, θmax=π, φmin=0, φmax=2π, meta...)
        nθ = round(Int, θmax / angres) + 1
        nφ = round(Int, φmax / angres) + 1
        dθ = θmax / nθ
        dφ = φmax / nφ
        θ = range(θmin, θmax, nθ)
        φ = range(φmin, φmax, nφ)
        new(PointsMonitor(name, origin; meta...), frame, radius, θ, φ)
    end
end
mutable struct ElectricalMonitor <: AbstractMonitor
    base_monitor::BaseMonitor
    voltage_line
    current_loop
    function ElectricalMonitor(name, voltage_line, current_loop; meta...)
        origin = (voltage_line[1] + voltage_line[2]) / 2
        # voltage_line=stack(collect.(voltage_line)) #
        # current_loop=stack(collect.(current_loop)) #
        new(BaseMonitor(name, origin, meta), voltage_line, current_loop)
    end
end

abstract type AbstractMonitorInstance end
@getr AbstractMonitorInstance

mutable struct PlaneMonitorInstance <: AbstractMonitorInstance
    name
    origin
    n
    dA
    A

    frame
    Is
    box_deltas
    local_plane_Is
    dr
    modes
    _modes
    λmodes
    _λmodes

    voltage_sample_Is
    current_sample_Iss
    voltage_dl
    current_dls
    λZ
end
isuniform(::PlaneMonitorInstance) = false

struct PointsMonitorInstance <: AbstractMonitorInstance
    # base_monitor::BaseMonitorInstance
    name
    origin
    points
    As
    Is
end
@functor PointsMonitorInstance (As,)


mutable struct SphereMonitorInstance <: AbstractMonitorInstance
    points_monitor_instance::PointsMonitorInstance
    n
    dA
    A

    θ
    φ
    frames
    invframes

    modes
    _modes
    λmodes
    _λmodes
end
isuniform(::SphereMonitorInstance) = true
Base.size(m::SphereMonitorInstance) = (length(m.θ), length(m.φ))
Base.length(m::SphereMonitorInstance) = length(m.θ) * length(m.φ)

@functor PlaneMonitorInstance (λmodes, _λmodes, dr, n, modes, _modes)
Base.ndims(m::PlaneMonitorInstance) = ndims(m.local_plane_Is)
area(m::PlaneMonitorInstance) = m.v
wavelengths(m::PlaneMonitorInstance) = keys(m.λmodes)
frequencies(m::PlaneMonitorInstance) = 1 ./ reverse(wavelengths(m))

@functor SphereMonitorInstance (points_monitor_instance, n, modes, _modes, λmodes,
    _λmodes, dA, points, frames, invframes, θ, φ)

mutable struct ElectricalMonitorInstance <: AbstractMonitorInstance
    voltage_sample_Is_monitor_instance::PointsMonitorInstance
    current_points_monitor_instance::PointsMonitorInstance
end
const MODES = ((Ex=1, Hy=1),
    (Ey=1, Hx=-1),
    (Ex=1, Hy=-1),
    (Ey=1, Hx=1),)

function _modes_from_n(n)
    # H, V, L, R polarizations
    (
        (Ex=1, Hy=n),
        (Ey=1, Hx=-n),
        (Ex=1, Ey=-im * 1, Hx=im * n, Hy=n),
        (Ex=1, Ey=im * 1, Hx=-im * n, Hy=n)),
    (
        (Ex=1, Hy=-n),
        (Ey=1, Hx=n),
        (Ex=1, Ey=-im * 1, Hx=-im * n, Hy=-n),
        (Ex=1, Ey=im * 1, Hx=im * n, Hy=-n)) ⊘ ((sqrt.(n),),)
end

function setup_line_integral(a, b, dr)
    T = eltype(a)
    l = b - a
    n = maximum(ceil.(Int, abs.(l) ./ dr))
    dl = l / n
    ps = [a + dl * (i - T(0.5)) for i = 1:n]
    # println("line integral with $n points")
    ps, dl
end
do_line_integral(u, ps, dl) =
    sum(ps) do p
        dl ⋅ map(_values(u)) do a
            a(p...)
        end
    end
function MonitorInstance(m::PlaneMonitor, λ, λs, g, ϵ; z=nothing, mode_solutions=nothing)
    @unpack F, bg = g

    i = findfirst(mode_solutions) do ms
        isnothing(ms.ports) || string(m.name) ∈ string.(ms.ports)
    end
    mode_solution = mode_solutions[i]
    @unpack λmodes, _λmodes, box_rulers, bbox, box_deltas, Is, dr, local_plane_Is, nmode, notPEC, dA, dr, plane_size = _get_stuff(m, λ, λs, ϵ, mode_solution, g; z)

    if isa(dA, AbstractArray)
        A = sum(dA)
    else
        A = dA * prod(size(local_plane_Is))
    end

    # heatmap(n) |> display
    # error()
    @unpack name, frame, start, origin, stop = m
    start, origin, stop = (start, origin, stop) / λ
    frame, start, origin, stop = FF.((frame, start, origin, stop))

    modes, _modes = _modes_from_n(nmode)
    modes = map(modes) do v
        vmap(v) do a
            a .* notPEC
        end
    end
    _modes = map(_modes) do v
        vmap(v) do a
            a .* notPEC
        end
    end

    @unpack voltage_line, current_loop = mode_solution
    voltage_sample_Is = voltage_dl = current_sample_Iss = current_dls = λZ = nothing
    sample_rulers = dr .* (:).(0, plane_size)
    _indexof(args...) = indexof.(args...)
    if !isnothing(voltage_line)
        v = _indexof.((sample_rulers,), voltage_line .- (start,)) - FF(0.5)
        voltage_sample_Is, v = setup_line_integral(v..., 1)
        voltage_dl = v .* dr
    end
    if !isnothing(current_loop)
        v = map(eachindex(current_loop)) do i
            p1 = current_loop[i]
            p2 = current_loop[mod1(i + 1, length(current_loop))]
            p1 = _indexof(sample_rulers, p1 - start) - FF(0.5)
            p2 = _indexof(sample_rulers, p2 - start) - FF(0.5)
            setup_line_integral(p1, p2, 1)
        end
        current_sample_Iss = first.(v)
        current_dls = [v .* dr for v = last.(v)]
    end
    if !isnothing(current_loop) && !isnothing(voltage_line)
        nλ = length(λs)
        λZ = zeros(complex(FF), nλ)
        for (i, d) = enumerate(λmodes[:, 1])
            V = do_line_integral([d(:Ex), d(:Ey)], voltage_sample_Is, voltage_dl)
            I = sum(zip(current_sample_Iss, current_dls)) do (pts, dl)
                do_line_integral([d(:Hx), d(:Hy)], pts, dl)
            end
            λZ[i] = V ./ I
        end
    end
    # @debug λZ
    r = PlaneMonitorInstance(name, origin, nmode, dA, A, frame, Is, box_deltas, local_plane_Is, dr, modes, _modes, λmodes, _λmodes,
        voltage_sample_Is, current_sample_Iss, voltage_dl, current_dls, λZ)
    r
end


function MonitorInstance(m::PointsMonitor, λ, points, rulers, offsets; kwargs...)
    @unpack name, origin = m
    origin /= λ
    v = kvmap(offsets) do k, v
        k => begin
            Is = map(points) do p
                indexof.(rulers, p,) - v
            end
            linix = LinearIndices(Tuple(length.(rulers) - 1))
            i_jv = reduce(vcat, map(enumerate(vec(Is))) do (i, _J)
                map(vec(Porcupine.nn(_J))) do (_J, w)
                    i, linix[_J...], w
                end
            end)
            _J = sort(getindex.(i_jv, 2) |> unique)
            d = Dict()
            for (j, _j) = enumerate(_J)
                d[_j] = j
            end
            A = sparse(getindex.(i_jv, 1), getindex.((d,), getindex.(i_jv, 2)), getindex.(i_jv, 3), length(points), length(_J))
            _J, A
        end
    end
    Is = vmap(x -> (x[1],), v)
    As = vmap(x -> x[2], v)
    PointsMonitorInstance(name, origin, points, As, Is)
end

function MonitorInstance(m::SphereMonitor, λ, λs, g, ϵ; mode_solutions=nothing, kwargs...)
    @unpack name, origin, radius, frame, θ, φ, meta = m
    @unpack F, bg, bbox, offsets, rulers = g
    origin, radius, frame, θ, φ, λ = FF.((origin, radius, frame, θ, φ, λ))
    points_monitor = PointsMonitor(name, origin; meta...)

    i = findfirst(mode_solutions) do ms
        isnothing(ms.ports) || string(m.name) ∈ string.(ms.ports)
    end
    mode_solution = isnothing(i) ? nothing : mode_solutions[i]

    dθ = step(θ)
    dφ = step(φ)
    R = radius / λ
    origin /= λ

    v = map(Base.product(θ, φ)) do (θ, φ)
        θ = max(FF(1f-3), θ)
        x = sin(θ) * cos(φ)
        y = sin(θ) * sin(φ)
        z = cos(θ)
        rhat = [x, y, z]

        r = R * rhat
        p = frame * r + origin |> collect

        θhat = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
        φhat = [-sin(φ), cos(φ), 0]


        _frame = [φhat -θhat rhat]
        _invframe = inv(_frame)

        p, _frame, _invframe
    end
    v = vec(v)
    points = getindex.(v, 1)
    frames = stack(getindex.(v, 2))
    invframes = stack(getindex.(v, 3))

    # n = sqrt.(samplemesh.((ϵ,), points; bg=bg(:ϵ)))
    n = 1
    modes, _modes = _modes_from_n(n)
    dA = R^2 * sin.(θ) * ones(FF, length(φ))' * dθ * dφ
    A = sum(dA)

    points_monitor_instance = MonitorInstance(points_monitor, λ, points, rulers, offsets)

    λmodes, _λmodes = _get_λsmodes(mode_solution, m, λ, λs, 1, dA)
    SphereMonitorInstance(points_monitor_instance, n, dA, A, θ, φ, frames, invframes, modes, _modes, λmodes, _λmodes)
    # error()
end
function MonitorInstance(m::ElectricalMonitor, g; kwargs...)
    @unpack voltage_line, current_loop, base_monitor = m
    @unpack F, bg, bbox, offsets, rulers = g
    @unpack name, origin = base_monitor
    voltage_line, current_loop, origin = FF.((voltage_line, current_loop, origin))

    points_monitor = PointsMonitor(name, origin; meta...)
    voltage_sample_Is_monitor_instance = MonitorInstance(points_monitor, λ, voltage_line, rulers, offsets)
    current_points_monitor_instance = MonitorInstance(points_monitor, λ, current_loop, rulers, offsets)

    ElectricalMonitorInstance(voltage_sample_Is_monitor_instance, current_points_monitor_instance)
    # error()
end

field(u::Map, k, m::AbstractMonitorInstance) = u(k)(m.Is[k]...)

function group_fields(u)
    E = u(r"E(.)")
    H = u(r"H(.)")
    (; E, H)
end

function energy(u, ϵ=1, Δ=1, Ibbox=nothing)
    if !isnothing(Ibbox)
        u, ϵ, Δ = bboxed.((u, ϵ, Δ), (Ibbox,))
    end
    mean(((ϵ .* .!(isPEC.(ϵ))) .* Porcupine.norm2(u.E) + Porcupine.norm2(u.H)) .* Δ)
end

_energy(u) =
    sum(u) do a
        mean(abs, a)
    end
