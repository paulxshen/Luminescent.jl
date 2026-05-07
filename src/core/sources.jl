# struct GaussianPulse
#     f
#     τ
#     d
#     b
# end
# duration(m::GaussianPulse) = m.d * m.τ
# fgp(t, f, τ, d, b) = (t <= d * τ) * (sqrt(2^-(2t / τ - d)^2) * cispi(2f * t) + b)
# function GaussianPulse(f::T, df) where {T}
#     f = T(f)
#     τ = T(0.44 / df)
#     d = T(2)
#     b = -mean(fgp.(range(0, d * τ, 10000), f, τ, d, 0))
#     GaussianPulse(f, τ, d, b)
# end
# (m::GaussianPulse)(t) = fgp(t, m.f, m.τ, m.d, m.b)

# struct SincPulse
#     f
#     f0
#     d
#     b
#     T
# end
# duration(m::SincPulse) = m.T
# fsp(t, f, f0, d, T, b) = (t <= T) * (sinc(2 * f0 * (t - d)) * cispi(2 * f * (t - d)) + b)
# function SincPulse(f::FF, df) where {FF}
#     f0 = FF(df / 2)
#     d = FF(1 / f0) / 2
#     T = 2d
#     b = -mean(fsp.(range(0, 2d, 10000), f, f0, d, T, 0))
#     SincPulse(f, f0, d, b, T)
# end
# (m::SincPulse)(t) = fsp(t, m.f, m.f0, m.d, m.T, m.b)

struct Pulse
    f
    T
    b
end
duration(m::Pulse) = m.T
_Pulse(t, f, T) = (t <= T) * cispi(2 * f * (t - T / 2)) * exp(-((t - T / 2) / (T / 4))^2)
function Pulse(f::FF, T::FF)
    b = -mean(_Pulse.(range(0, T, 1000), f, T))
    # b = 0
    Pulse(f, T, b)
end
(m::Pulse)(t) = _Pulse(t, m.f, m.T) + m.b

struct Wave
    f
end
duration(m::Wave) = Inf
(m::Wave)(t) = cispi(2m.f * t)

function completexyz(mode, N)
    # length(mode) == N && return mode
    OrderedDict([
        begin
            i = findfirst(keys(mode)) do k
                endswith(string(k), s)
            end
            f = first(string(first(keys(mode))))
            # @show keys(mode), s, i
            if isnothing(i)
                (Symbol("$f$s") => 0)
            else
                (keys(mode)[i] => mode[keys(mode)[i]])
            end
        end for s = "xyz"[1:N]
    ])
end

# """
# """
mutable struct Source
    wavelength
    bandwidth
    duration
    modenums

    origin
    start
    stop
    frame
    # dimsperm

    name
    meta
end
@getr Source
function Source(wavelength, bandwidth, duration, origin, start, stop, frame; modenums=[0], name, meta...)
    Source(wavelength, bandwidth, duration, modenums, origin, start, stop, frame, name, meta)
end

Base.ndims(m::Source) = length(m.origin)


struct SourceInstance
    name
    origin
    sigmodes
    meta
end
@getr SourceInstance
@functor SourceInstance (sigmodes,)
Base.ndims(m::SourceInstance) = length(m.origin)
duration(m::SourceInstance) = maximum(duration.(first.(m.sigmodes)))

function SourceInstance(s::Source, λ, λs, g, ϵ; z=nothing, mode_solutions=nothing)
    @unpack name, frame, origin, modenums, wavelength, bandwidth, duration, meta = s
    @unpack F, offsets, sz, dt, weights = g

    name = "source_$name"
    wavelength, bandwidth, duration, origin = FF.((wavelength, bandwidth, duration, origin))
    wl = wavelength / λ
    origin /= λ
    N = ndims(s)


    i = findfirst(mode_solutions) do ms
        isnothing(ms.ports) || string(s.name) ∈ string.(ms.ports)
    end
    mode_solution = mode_solutions[i]
    @unpack λmodes, _λmodes, box_size, bbox, box_deltas, box_weights, Is, local_plane_Is, dA =
        _get_stuff(s, λ, λs, ϵ, mode_solution, g; z)
    modes = λmodes[findmin(abs, λs - wl)[2], :]

    if isnothing(bandwidth)
        f = λ / (wavelength)
        f = Pulse(f, duration)
    else
        f1 = λ / (wavelength + bandwidth / 2)
        f2 = λ / (wavelength - bandwidth / 2)
        df = f2 - f1
        f = (f1 + f2) / 2
        if f == 0
            Wave(f)
        else
            # @show f, df, dt
            f = Pulse(f, dt * round(1 / FF(df) / dt) |> FF)
        end
    end

    global sigmodes = tuple.((f,), modes[modenums+1])
    # f = SincPulse(1 / wl, 1.5df)
    # y = abs(fft(f.(0:0.01:50)))
    # x = range(0, 100, length(y))
    # n = round(Int, length(y) / 20)
    # plot(x[1:n], y[1:n]) |> display
    # error("Not implemented yet")

    sigmodes = [
        begin
            mode = vmap(mode) do v
                # Symbol(k) == :Jy && (heatmap(mode(k) |> real) |> display) && error("stop") 
                a = zeros(complex(FF), box_size)
                for (I, w, v, dA) = tuple.(local_plane_Is, box_weights, v, dA)
                    I, w, v, dA = (I, w, v, dA) |> FF
                    place!(a, v / w * dA, I)
                end
                # Symbol(k) == :Jy && (heatmap(a[1, :, :] |> real) |> display) && error("stop")
                a
            end

            mode = completexyz(mode, N)

            mode = todvoa(mode)
            mode = vmap(mode) do v
                frame * v[1:N]
            end
            mode = todxyz(mode)
            mode = namedtuple([
                k => begin
                    v = mode(k)
                    @assert !any(isnan, v)

                    if all(iszero, v)
                        a = 0
                    else
                        a = zeros(complex(FF), sz)
                        place!(a, v, first.(Is[k]))
                    end
                    a |> FF
                end for k = sort(keys(mode))])
            @debug keys(mode)
            (f, mode)
        end for (f, mode) = sigmodes
    ]
    SourceInstance(name, origin, sigmodes, meta)
end

# Complex
function (s::SourceInstance)(t::FF)
    namedtuple([
        k => sum(s.sigmodes) do (f, _mode)
            t > duration(f) && return 0
            real(complex(FF)(f(t)) .* _mode[k])
        end for k = keys(s.sigmodes[1][2])])
end

function E2J(d)
    k0 = filter(k -> string(k)[1] == 'E', keys(d))
    k = [Symbol("J" * string(k)[2:end]) for k in k0]
    v = [d[k] for k in k0]
    dict(Pair.(k, v))
end
function EH2JM(d::T) where {T}
    dict(Pair.(replace(keys(d), :Ex => :Jx, :Ey => :Jy, :Ez => :Jz, :Hx => :Mx, :Hy => :My, :Hz => :Mz), values(d)))
end

function _get_λsmodes(mode_solution, sm, λ, λs, ϵmode)
    @unpack dA, plane_size = mode_solution
    isnothing(mode_solution) && return nothing, nothing
    # c = isa(dA, Number) ? 1 / sqrt(length(ϵmode)) : 1
    @unpack modes, wavelengths = mode_solution
    # @debug wavelengths
    a, b = extrema(wavelengths)
    # @debug size(modes)

    nλ = length(wavelengths)
    x = (2wavelengths - a - b) ./ (b - a)
    M = reduce(hcat, [x .^ i for i = 0:nλ-1]) |> inv
    C = map(eachcol(modes)) do modes
        map(eachrow(M)) do v
            sum(v .* modes)
        end
    end |> stack

    λmodes = permutedims(stack(map(λs) do w
        # x = (2w - a - b) / (b - a)
        # map(eachcol(C)) do v
        #     normalize_mode(sum(v .* (x .^ (0:nλ-1))) |> FF, dA) |> FF
        # end
        # i = argmin(1:nλ) do i
        #     abs(wavelengths[i] - w)
        # end
        # modes[i, :]

        i = searchsortedfirst(wavelengths, w)
        i == 1 && return modes[1, :]
        i > length(wavelengths) && return modes[end, :]
        a, b = (wavelengths[i-1], wavelengths[i])
        wa, wb = (b - w) / (b - a), (w - a) / (b - a)
        modes[i-1, :] * wa + modes[i, :] * wb
    end))

    if sm isa PlaneMonitor || sm isa Source
        λmodes = map(λmodes) do mode
            normalize_mode(mode, dA)
        end
    end
    # global b = λmodes
    if isa(sm, AbstractMonitor)
        _λmodes = mirror_mode.(λmodes)
    else
        λmodes = addJ.(λmodes, (ϵmode,))
        _λmodes = nothing
    end

    λmodes, _λmodes
end
function _get_stuff(sm, λ, λs, ϵ, mode_solution, g; z)
    @unpack origin, start, stop, frame, name = sm
    @unpack plane_size, dr, dA = mode_solution
    @unpack F, rulers, offsets, weights, bg, searchers, funnel = g
    bg = bg(:ϵ)
    λ, origin, start, stop, frame = FF.([λ, origin, start, stop, frame])

    origin /= λ
    start /= λ
    stop /= λ
    P = frame[:, 1:end-1]
    center = origin + P * (start + stop) / 2

    if isa(sm, Source)
        ks = [k for k = keys(offsets) if string(k)[1] ∈ ('J')]
    elseif isa(sm, AbstractMonitor)
        ks = [k for k = keys(offsets) if string(k)[1] ∈ ('E', 'H')]
    end

    s = Int.(sign.(P * (stop - start)))
    start = indexof.(rulers, origin + P * start)
    stop = indexof.(rulers, origin + P * stop)

    s += s .== 0
    start = s .* (floor.(Int, s .* start))
    stop = s .* (ceil.(Int, s .* stop))
    box_size = Tuple(max.(1, abs.(stop - start)))

    signed_plane_rulers = getindex.(rulers, map(start, stop) do a, b
        a == b ? (a:b) : a:sign(b - a):b
    end)
    signed_plane_rulers -= first.(signed_plane_rulers)

    start, stop = min.(start, stop), max.(start, stop)
    bbox = [getindex.(rulers, start) getindex.(rulers, stop)]
    box_rulers = getindex.(rulers, range.(start, stop))
    box_deltas = diff.(box_rulers)

    Is = dict([k => begin
        a = start - offsets[k]
        a = map(a) do r
            if r < 2
                1
            else
                r
            end
        end
        range.(a, a + box_size - 1, box_size)
    end for k = ks])

    signed_plane_start = center - P * (FF(collect((plane_size - 1) / 2)) .* dr)

    global_plane_points = map(CartesianIndices(Tuple(plane_size))) do I
        P * (collect(Tuple(I) - 1) .* dr) + signed_plane_start |> FF
    end
    global_plane_Is = map(global_plane_points) do p
        indexof.(rulers, p) - 0.5 |> FF
    end
    box_weights = (splat(weights)).(global_plane_Is)
    local_plane_Is = map(global_plane_Is) do v
        v = v - start + 1
        v = max.(v, 1)
        min.(v, box_size) |> Tuple
    end

    tol = mean(mean.(box_deltas)) / 100
    ϵmode = samplemesh.(map(v -> length(v) == 3 ? v : SVector{3}(v..., z), global_plane_points), (funnel,), (searchers,), (ϵ,), (tol,)) .|> FF
    # heatmap(ϵmode) |> display
    # error("stop")

    λmodes, _λmodes = _get_λsmodes(mode_solution, sm, λ, λs, ϵmode)

    notPEC = .!(isPEC(ϵmode))
    nmode = sqrt.(ϵmode)

    r = (; λmodes, _λmodes, box_size, box_rulers, bbox, box_deltas, Is, box_weights, local_plane_Is, plane_size, ϵmode, nmode, notPEC, dA, dr)
end