function resize(a, sz)
    size(a) != Tuple(sz) && @warn "array size $(size(a)) not same as new size $sz. array will be interpolated"
    imresize(a, sz, method=ImageTransformations.Lanczos4OpenCV())

end

function _aug(mode, N)
    length(mode) == N && return mode
    OrderedDict([
        begin
            i = findfirst(keys(mode)) do k
                endswith(string(k), s)
            end
            # @show keys(mode), s, i
            isnothing(i) ? (Symbol("D$s") => 0) : (keys(mode)[i] => mode[keys(mode)[i]])
        end for s = "xyz"[1:N]
    ])
end

"""
"""
struct Source
    specs
    center
    lb
    ub
    # normal
    # tangent
    # zaxis
    # xaxis
    dimsperm
    N
    center3
    lb3
    ub3
    tags
    # function Source(sigmodes, center::Base.AbstractVecOrTuple, normal, tangent, lb, ub, ; tags...)
    #     new(f, center, lb, ub, normal, tangent, E2J(mode), tags)
    # end
    # function Source(sigmodes, center::Base.AbstractVecOrTuple, L; tags...)
    #     new(sigmodes, center, -L / 2, L / 2, getdimsperm(L), tags)
    # end

end
Source(args...; λmodenums=nothing, λmodes=nothing, tags...) = Source((; λmodenums, λmodes), args..., tags)
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
    o
    center
    dimsperm
    ϵeff
    tags
end
@functor SourceInstance (sigmodes,)

function SourceInstance(s::PlaneWave, g)
    @unpack L = g
    @unpack dims, sigmodes, tags = s
    SourceInstance(Source(sigmodes, L / 2, -L / 2, L / 2, getdimsperm(dims), tags), g)
end

function SourceInstance(s::Source, g, ϵ, temp, mode_solutions=nothing)
    @unpack center, lb, ub, tags, dimsperm, specs, N, center3, lb3, ub3 = s
    @unpack F, deltas, deltas3, field_sizes, field_lims, mode_spacing, dl = g
    C = complex(F)
    ϵeff = nothing
    dx = deltas[1][1]
    @unpack λmodenums, λmodes = specs
    if !isnothing(λmodenums)
        start = int((center3 + lb3) / dl)
        stop = int((center3 + ub3) / dl)
        start, stop = min.(start, stop), max.(start, stop)

        sel = abs.(stop - start) .>= 0.001
        stop[!sel] .= start[!sel]
        start += 0.5sel
        stop -= 0.5sel

        start += 0.5
        stop += 0.5
        len = int(stop - start + 1)
        ϵmode = getindexf(ϵ, range.(start, stop, len)...)
        ϵmode = permutedims(ϵmode, dimsperm, 2)

        # global _a = ϵmode, dimsperm
        λmodes = OrderedDict([λ => begin
            modes = solvemodes(ϵmode, dl, λ, maximum(mns) + 1, mode_spacing, temp; mode_solutions)[mns+1]
            if isnothing(ϵeff)
                Ex = modes[1].Ex
                v = real(sum(downsample(ϵmode, mode_spacing) .* Ex, dims=2) ./ sum(Ex, dims=2))
                ϵeff = maximum(ϵmode) => v[round(length(v) / 2)]
            end
            modes
        end for (λ, mns) = pairs(λmodenums)])
        if N == 2
            λmodes = kmap(v -> collapse_mode.(v), λmodes)
        end
    end

    start = v2i(center + lb - g.lb, deltas)
    stop = v2i(center + ub - g.lb, deltas)
    start, stop = min.(start, stop), max.(start, stop)

    sel = abs.(stop - start) .>= 0.001
    stop[!sel] .= start[!sel]
    start += 0.5sel
    stop -= 0.5sel

    o = NamedTuple([k => F.(start - fl[:, 1] + 1) for (k, fl) = pairs(field_lims)])
    λmodes = fmap(F, λmodes)
    sigmodes = reduce(vcat, [collect(zip(fill(λ, length(modes)), modes)) for (λ, modes) = (pairs(λmodes))])
    sigmodes = [
        begin
            _f = if isa(sig, Number)
                t -> cispi(2t / sig)
            else
                sig
            end
            f = x -> convert(C, _f(x))

            mode = NamedTuple([k => v for (k, v) = pairs(mode) if startswith(string(k), "D")])

            mode = permutexyz(mode, invperm(dimsperm), N)
            mode = _aug(mode, N)
            ks = sort([k for k = keys(mode) if string(k)[end] in "xyz"[1:N]])
            _mode = namedtuple([k => begin
                b = mode(k)
                k = Symbol("E$(string(k)[end])")
                a = zeros(C, field_sizes[k])
                if b == 0
                    0
                else
                    # global aaaa = a, mode, o, k
                    setindexf!(a, b, range.(o[k], o[k] + size(b) - 1)...)
                    a
                end
            end for k = ks])
            (f, _mode)
        end for (sig, mode) = sigmodes
    ]
    # error("stop")
    _center = round(v2i(center - g.lb, deltas) + 0.5)
    # @show center, g.lb, _center

    SourceInstance(sigmodes, o, _center, dimsperm, ϵeff, tags)
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
# function EH2JM(d::NamedTuple)
#     NamedTuple((d) |> pairs |> OrderedDict |> EH2JM |> pairs)
# end


# function SourceInstance(s::PlaneWave, deltas, field_sizes, common_left_pad_amount, fl, sz0; F=Float32)
#     @unpack sigmodes, dims,  tags = s
#     _F(x::Real) = F(x)
#     _F(x::Complex) = ComplexF32(x)
#     f = _F ∘ f
#     d = length(common_left_pad_amount)
#     g = Dict([k => _F.(mode[k]) * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / deltas for k = keys(mode)])
#     o = NamedTuple([k =>
#         1 .+ fl[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? field_sizes[k][i] - 1 : 0 for i = 1:d])
#                     for k = keys(fl)])
#     _mode = Dict([k => place(zeros(F, field_sizes[k]), o[k], mode[k],) for k = keys(mode)])
#     c = first(values(field_sizes)) .÷ 2
#     SourceInstance(f, g, _mode, o, c,  tags)
# end

# function SourceInstance(s::GaussianBeam, deltas, field_sizes, fl, stop; F=Float32)
#     _F(x::Real) = F(x)
#     _F(x::Complex) = ComplexF32(x)
#     f = _F ∘ f
#     @unpack f, σ, mode, c, dims = s
#     n = round(Int, 2σ / deltas)
#     r = n * deltas
#     r = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(c)]
#     g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(r...)] / deltas
#     fl = fl .- 1 .+ index(c, deltas) .- round.(Int, (size(g) .- 1) ./ 2)
#     _mode = place(zeros(F, sz), g, fl)
#     SourceInstance(f, g, _mode, fl, c, tags)
# end



# J = OrderedDict()
# for k = (:Jx, :Jy, :Jz)
#     if k in keys(s.mode)
#         J[k] = ComplexF32.(s.mode[k])
#     else
#         J[k] = zeros(ComplexF32, size(s.mode(1)))
#     end
# end
# J = values(J)

# L = ub .- lb
# d = length(lb)
# D = length(center) # 2D or 3D
# if D == 2
#     zaxis = [normal..., 0]
#     yaxis = [0, 0, 1]
#     # xaxis = cross(yaxis, zaxis)
#     xaxis = cross(yaxis, zaxis)
# else
#     zaxis = convert.(F, normal)
#     xaxis = convert.(F, tangent)
#     yaxis = cross(zaxis, xaxis)
# end

# frame = [xaxis, yaxis, zaxis]
# J = reframe(frame, J)
# if D == 2
#     J = [J[:, :, 1] for J in J]
# end
# lb = sum(lb .* frame[1:d])[1:D]
# ub = sum(ub .* frame[1:d])[1:D]
# L = ub - lb

# Jx, Jy, Jz = J
# if D == 2
#     mode = (; Jx, Jy)
# else
#     mode = (; Jx, Jy, Jz)
# end
# v = zip(lb, ub)
# lb = minimum.(v)
# ub = maximum.(v)
# n = findfirst(abs.(zaxis) .> 0.001)
# mode = ignore_derivatives() do
#     mode / deltas[n][1]#[findfirst(>(center[n]), cumsum(deltas[n]))]
# end
