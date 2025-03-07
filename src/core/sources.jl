
function completexyz(mode, N)
    # length(mode) == N && return mode
    OrderedDict([
        begin
            i = findfirst(keys(mode)) do k
                endswith(string(k), s)
            end
            F = first(string(first(keys(mode))))
            # @show keys(mode), s, i
            if isnothing(i)
                (Symbol("$F$s") => 0)
            else
                (keys(mode)[i] => mode[keys(mode)[i]])
            end
        end for s = "xyz"[1:N]
    ])
end

"""
"""
mutable struct Source
    λmodenums
    λsmode
    λmodes

    center
    dimensions
    frame
    # dimsperm

    tags
end

Source(center, dimensions, frame; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...) =
    Source(λmodenums, λsmode, λmodes, center, dimensions, frame, tags)

Base.ndims(m::Source) = length(m.center)
# isortho(m::Source) = !isnothing(m.dimsperm)

"""
    function PlaneWave(f, dims; mode...)

Constructs plane wave source

Args
- f: time function
- dims: eg -1 for wave coming from -x face
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct PlaneWave
    λmodes
    dims
    tags
    function PlaneWave(sigmodes, dims; tags...)
        new(sigmodes, dims, tags)
    end
end
@functor PlaneWave
"""
    function GaussianBeam(f, σ, c, dims; mode...)

Constructs gaussian beam source

Args
- f: time function
- σ: std dev length
- dims: eg 1 for x direction
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct GaussianBeam
    f
    σ
    mode
    c
    dims
    function GaussianBeam(f, σ, c, dims; mode...)
        new(f, σ, mode, c, dims)
    end
end
@functor GaussianBeam

"""
    function Source(f, c, lb, ub, tags=""; mode...)
    function Source(f, c, dimensions, tags=""; mode...)
        
Constructs custom  source. Can be used to specify uniform or modal sources

Args
- f: time function
- c: g.lb or center of source
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c
- dimensions: source dimensions in [wavelengths]
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""

# mode(m::Source) = m.mode
# Base.string(m::Union{Source,Source}) =
#     """
#     $(m.tags): $(count((m.ub.-m.lb).!=0))-dimensional source, centered at $(m.center|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center,"""
#  exciting $(join(keys(m.mode),", "))"""


struct SourceInstance
    sigmodes
    # o
    center
    # dimsperm
    tags
end
@functor SourceInstance (sigmodes,)
Base.ndims(m::SourceInstance) = length(m.center)

function SourceInstance(s::PlaneWave, g)
    @unpack dimensions = g
    @unpack dims, sigmodes, tags = s
    SourceInstance(Source(sigmodes, dimensions / 2, -dimensions / 2, dimensions / 2, getdimsperm(dims), tags), g)
end

function SourceInstance(s::Source, g, ϵ, TEMP, mode_solutions=nothing)
    @unpack dimsperm, tags = s
    @unpack F, field_sizes = g
    C = complex(F)
    N = ndims(s)
    ϵeff = nothing

    λmodes, _λmodes, plane_rulers, bbox, plane_deltas, I, plane_points, labelpos = _get_λmodes(s, ϵ, TEMP, mode_solutions, g)
    # global _a = λmodes

    λs = @ignore_derivatives Array(keys(λmodes))
    modess = values(λmodes)
    if all(x -> x === (modess[1]), modess)
        iss = [eachindex(λs)]
        println("all modes are the same")
    else
        iss = cluster(λs)
    end
    sigmodes = map(iss) do is
        f = t -> sum(getindex.((Array(λs),), Array(is))) do λ
            cispi(2t / λ) |> C
        end
        _modess = getindex.((modess,), is)
        modes = _modess[round(length(_modess) / 2 + 0.1)]
        f, modes
    end
    global _a1 = sigmodes

    sigmodes = reduce(vcat, [collect(zip(fill(f, length(modes)), modes)) for (f, modes) = sigmodes])

    sigmodes = [
        begin
            _f = if isa(sig, Number)
                t -> cispi(2t / sig)
            else
                sig
            end
            f = x -> convert(C, _f(x))

            mode = OrderedDict([
                begin
                    a = zeros(F, length.(plane_deltas))
                    for (p, v) = zip(plane_points, v)
                        I = indexof.(plane_rulers, p) - F(0.5)
                        place!(a, v, I)
                    end
                    k => a
                end for (k, v) = pairs(mode)])
            mode = completexyz(mode, N)
            # ks = sort([k for k = keys(mode) if string(k)[end] in "xyz"[1:N]])
            # println(ks)
            mode = namedtuple([
                begin
                    v = mode(k)
                    @assert !any(isnan, v)

                    if all(iszero, v)
                        v = 0
                    else
                        a = zeros(C, field_sizes[k])
                        place!(a, v, first.(I[k]))
                        v = a
                    end
                    @assert !any(isnan, v)
                    k => v
                end for k = sort(keys(mode))])

            (f, mode)
        end for (sig, mode) = sigmodes]
    # @show center, g.lb, labelpos

    SourceInstance(sigmodes, labelpos, tags)
end

# Complex
function (s::SourceInstance)(t::Real)
    namedtuple([
        k => sum(s.sigmodes) do (f, _mode)
            real(f(t) .* _mode[k])
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

function _get_λmodes(s, ϵ, TEMP, mode_solutions, g)
    @unpack center, dimensions, tags, frame, λmodenums, λmodes, λsmode = s
    @unpack F, deltas, deltas3, field_sizes, field_lims, mode_spacing, spacings, dl, padamts = g
    N = ndims(s)

    C = complex(g.F)

    if isa(s, Source)
        ks = [k for k = keys(field_lims) if string(k)[1] ∈ ('J')]
    elseif isa(s, Monitor)
        ks = [k for k = keys(field_lims) if string(k)[1] ∈ ('E', 'H')]
    end

    start = floor.(Int, indexof.(rulers, center - dimensions / 2))
    stop = ceil.(Int, indexof.(rulers, center + dimensions / 2))
    bbox = [getindex.(rulers, start) getindex.(rulers, stop)]

    plane_rulers = getindex.(rulers, range.(start, stop))
    plane_deltas = diff.(plane_rulers)

    len = stop - start
    I = dict([k => begin
        l, r = eachcol(offsets[k])
        range.(start - l, stop - r - 1, len)
    end for k = ks])
    labelpos = round(Int, indexof.(rulers, center))

    # 
    if !isnothing(λmodenums)
        sz = round(dimensions / dx)
        P = frame[1:end-1, 1:end-1]
        plane_start = center - P * collect((sz - 1) / 2) - bbox[:, 1]
        plane_points = map(CartesianIndices(Tuple(sz))) do I
            P * collect(Tuple(I) - 1) + plane_start
        end
        plane_Is = map(plane_points) do p
            indexof.(plane_rulers, p)
        end

        ϵmode = samplemesh(ϵ, plane_points .+ bbox[:, 1])
        λmodes = OrderedDict([λ => begin
            modes = solvemodes(ϵmode, dx, λ, maximum(mns) + 1, TEMP; mode_solutions)[mns+1]
            map(modes) do mode
                if isa(s, Monitor)
                    ks = filter(k -> string(k)[1] ∈ ('E', 'H'), keys(mode))
                else
                    ks = filter(k -> string(k)[1] ∈ ('J'), keys(mode))
                end
                namedtuple(ks .=> mode.(ks))
            end
        end for (λ, mns) = pairs(λmodenums)])

        if isa(s, Monitor)
            λmodes = kmap(λmodes) do modes
                normalize_mode.(modes, (plane_deltas,))
            end
            _λmodes = kmap(λmodes) do modes
                mirror_mode.(modes)
            end
        else
            _λmodes = nothing
        end
    elseif !isnothing(λsmode)
        # λs, mode = λsmode
        # λs = F.(λs)
        # mode = kmap(Symbol, identity, mode)
        # mode = OrderedDict([
        #     begin
        #         v = mode[k]
        #         if isa(v, Number)
        #             v = v * block
        #         end
        #         k => v
        #     end for k = keys(mode) if k ∈ ks])

        # if isa(s, Monitor)
        #     mode = normalize_mode(mode, md)
        #     _mode = mirror_mode(mode; flip=false)
        # end

        # λmodes = OrderedDict([λ => [mode] for λ = λs])
        # if isa(s, Monitor)
        #     _λmodes = OrderedDict([λ => [_mode] for λ = λs])
        # else
        #     _λmodes = nothing
        # end
    end

    λmodes = fmap(F, λmodes)
    λmodes = sort(λmodes, by=kv -> kv[1])

    if !isnothing(_λmodes)

        _λmodes = fmap(F, _λmodes)
        _λmodes = sort(_λmodes, by=kv -> kv[1])
    end
    # m = kmap(x -> C.(x), m)

    (; λmodes, _λmodes, plane_rulers, bbox, plane_deltas, I, plane_points, plane_Is, labelpos)
end