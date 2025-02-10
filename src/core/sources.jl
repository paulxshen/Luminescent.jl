
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
    L
    center3
    L3
    mask
    # normal
    # tangent
    # zaxis
    # xaxis
    dimsperm
    frame

    approx_2D_mode
    tags
    # function Source(sigmodes, center::Base.AbstractVecOrTuple, normal, tangent, lb, ub, ; tags...)
    #     new(f, center, lb, ub, normal, tangent, E2J(mode), tags)
    # end
    # function Source(sigmodes, center::Base.AbstractVecOrTuple, L; tags...)
    #     new(sigmodes, center, -L / 2, L / 2, getdimsperm(L), tags)
    # end

end

Source(center, L, dimsperm, center3=center, L3=L, approx_2D_mode=nothing; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...) =
    Source(λmodenums, λsmode, λmodes, center, L, center3, L3, nothing, dimsperm, nothing, approx_2D_mode, tags)

Source(mask, frame; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...) =
    Source(λmodenums, λsmode, λmodes, nothing, nothing, nothing, nothing, mask, nothing, frame, nothing, tags)

function Base.ndims(m::Source)
    v = m.center
    if !isnothing(v)
        length(v)
    else
        ndims(m.mask)
    end
end

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
    function Source(f, c, L, tags=""; mode...)
        
Constructs custom  source. Can be used to specify uniform or modal sources

Args
- f: time function
- c: g.lb or center of source
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c
- L: source dimensions in [wavelengths]
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
    @unpack L = g
    @unpack dims, sigmodes, tags = s
    SourceInstance(Source(sigmodes, L / 2, -L / 2, L / 2, getdimsperm(dims), tags), g)
end

function SourceInstance(s::Source, g, ϵ, TEMP, mode_solutions=nothing)
    @unpack dimsperm, tags = s
    @unpack F, field_sizes = g
    C = complex(F)
    N = ndims(s)
    ϵeff = nothing

    λmodes, _, inds, labelpos = _get_λmodes(s, ϵ, TEMP, mode_solutions, g)
    global _a = λmodes

    λs = @ignore_derivatives Array(keys(λmodes))
    modess = values(λmodes)
    if !isnothing(s.frame)
        #  all(x -> x === (modess[1]), modess) 
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

            mode = globalframe(mode, s)
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
                        I = inds[k]
                        # global aaaa = a, b, mode, k, I
                        setindexf!(a, v, I...;)# approx=true)
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
    @unpack center, L, tags, dimsperm, frame, center3, L3, approx_2D_mode, λmodenums, λmodes, λsmode = s
    @unpack F, deltas, deltas3, field_sizes, field_lims, mode_spacing, spacings, dl, padamts = g
    N = ndims(s)

    C = complex(g.F)
    md = first.(g.mode_deltas)

    labelpos = zeros(F, N)
    if isa(s, Source)
        ks = [k for k = keys(field_lims) if string(k)[1] ∈ ('J')]
    elseif isa(s, Monitor)
        ks = [k for k = keys(field_lims) if string(k)[1] ∈ ('E', 'H')]
    end

    dx = deltas[1][1]

    start = v2i(center - L / 2 - g.lb, deltas)
    stop = v2i(center + L / 2 - g.lb, deltas)
    start, stop = min.(start, stop), max.(start, stop)

    sel = abs.(start - stop) .> 1e-3
    start += 0.5sel
    stop -= 0.5sel
    len = Tuple(int(stop - start + 1))

    start = F(start)
    stop = F(stop)

    inds = dict([k => begin
        lr = field_lims[k]
        range.(start, stop, len) - lr[:, 1] + 1
    end for k = ks])
    labelpos = round(v2i(center - g.lb, deltas) + 0.5)

    # 
    if !isnothing(λmodenums)

        start3 = round((center3 - L3 / 2) / dl + 0.001)
        stop3 = round((center3 + L3 / 2) / dl + 0.001)
        start3, stop3 = min.(start3, stop3), max.(start3, stop3)

        sel3 = start3 .!= stop3
        start3 += 0.5sel3
        stop3 -= 0.5sel3

        ratio = int(dx / dl)
        stop3[1:N] = sel3[1:N] .* (ratio * len[1:N] - 1) + start3[1:N]
        len3 = int(stop3 - start3 + 1)

        start3 += 0.5
        stop3 += 0.5

        start3 = F(start3)
        stop3 = F(stop3)

        ϵmode = getindexf(ϵ, range.(start3, stop3, len3)...;)
        # global _a = ϵmode, start3, stop3, len3, dimsperm
        ϵmode = permutedims(ϵmode, dimsperm, 2)

        # global _a = ϵmode, dimsperm
        λmodes = OrderedDict([λ => begin
            modes = solvemodes(ϵmode, dl, λ, maximum(mns) + 1, mode_spacing, TEMP; mode_solutions)[mns+1]
            # if isnothing(ϵeff)
            #     @unpack Ex, Ey = modes[1]
            #     E = sqrt.(abs2.(Ex) + abs2.(Ey))
            #     ϵ = downsample(ϵmode, mode_spacing)
            #     J = sum(ϵ .* E, dims=2)
            #     E = sum(E, dims=2)
            #     ϵ = J ./ E
            #     ϵmin, ϵmax = extrema(ϵmode)
            #     # v = filter(v) do x
            #     #     ϵmin <= x <= ϵmax
            #     # end

            #     # ϵeff = maximum(ϵmode) => ϵ[argmax(E)]
            #     # ϵeff = maximum(ϵmode) => 0.8ϵmax + 0.2ϵmin
            #     ϵeff = maximum(ϵmode) => ϵmax
            # end
            map(modes) do mode
                if N == 2
                    mode = collapse_mode(mode, approx_2D_mode)
                end
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
                normalize_mode.(modes, (md,))
            end
            _λmodes = kmap(λmodes) do modes
                mirror_mode.(modes)
            end
        else
            _λmodes = nothing
        end
    elseif !isnothing(λsmode)
        λs, mode = λsmode
        λs = F.(λs)
        mode = kmap(Symbol, identity, mode)
        mode = OrderedDict([
            begin
                v = mode[k]
                if isa(v, Number)
                    v = v * block
                end
                k => v
            end for k = keys(mode) if k ∈ ks])

        if isa(s, Monitor)
            mode = normalize_mode(mode, md)
            _mode = mirror_mode(mode; flip=false)
        end

        λmodes = OrderedDict([λ => [mode] for λ = λs])
        if isa(s, Monitor)
            _λmodes = OrderedDict([λ => [_mode] for λ = λs])
        else
            _λmodes = nothing
        end
    end

    λmodes = fmap(F, λmodes)
    λmodes = sort(λmodes, by=kv -> kv[1])

    if !isnothing(_λmodes)

        _λmodes = fmap(F, _λmodes)
        _λmodes = sort(_λmodes, by=kv -> kv[1])
    end
    # m = kmap(x -> C.(x), m)

    λmodes, _λmodes, inds, labelpos
end