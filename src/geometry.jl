struct SubpixelAveraging
    v
end

function _apply_subpixel_averaging(lr, a::AbstractArray)
    N = ndims(a)
    lr = cpu(lr)
    l, r = eachcol(lr)
    a = pad(a, :replicate, 1)
    _l = l + 1.5
    _r = _l + r - l
    getindexf(a, range.(_l, _r + 0.001)...)
end
#     # for (i, (do_smooth, pad_left, pad_right)) in enumerate(s.v)
#     for (i, (l, r)) = enumerate(eachcol(s.v))
#         select = i .== (1:d) |> Tuple
#         if l == -1
#             a = pad(a, :replicate, select, 0)
#         end
#         if r == 1
#             a = pad(a, :replicate, 0, select)
#         end
#         if l == -1 == r
#             a = (2selectdim(a, i, 2:(size(a, i))) - diff(a, dims=i)) / 2
#         elseif l == 1 == r
#             a = (2selectdim(a, i, 1:(size(a, i)-1)) + diff(a, dims=i)) / 2
#         end
#     end
#     a
# end


function apply_subpixel_averaging(sas, gs)
    namedtuple([k => _apply_subpixel_averaging(sas[k], values(gs)[findfirst(keys(gs)) do gk
        startswith(string(k), string(gk))
    end]) for k = keys(sas)])
end
function apply_tensor_smoothing(geomlims, A)
    (; invϵ=stack([_apply_subpixel_averaging.((gl,), c) for (c, gl) = zip(eachcol(A), geomlims(r"ϵ.*"))]),)
end

function _apply_geometry_padding(p::AbstractVector{<:OutPad}, a)
    for p = p
        @unpack l, r, b = p
        a = pad(a, b, l, r)
    end
    a
end

function apply_geometry_padding(gps, gs)
    namedtuple([k => k in keys(gps) ? _apply_geometry_padding(gps[k], gs[k]) : gs[k] for k = keys(gs)])
end

function _apply_field_padding(p::AbstractVector{<:InPad}, a::AbstractArray; nonzero_only=false)
    if nonzero_only
        p = filter(p) do p
            p.b != 0
        end
    end
    isempty(p) && return a

    if length(p) == 1 && !isnothing(p[1].m)
        return a .* p[1].m
    end

    a_ = Buffer(a)
    # a_ .= a
    a_[axes(a)...] = a
    for p = p
        @unpack l, r, b, m = p
        if isnothing(m)
            a_ = pad!(a_, b, l, r;)
        else
            a_[axes(a)...] = a .* m
        end
    end

    copy(a_)
end

function apply_field_padding(fps, fs; kw...)
    # namedtuple([k => apply(fps[k], fs[k]; kw...) for k = keys(fs)])
    dict([k => k in keys(fps) ? _apply_field_padding(fps[k], fs[k]; kw...) : fs[k]
          for k = keys(fs)])
end

function tensorinv(a, ratio)
    N = ndims(a)
    T = eltype(a)
    inva = downsample(a, ratio) do a
        n = zeros(T, N)
        if maximum(a) != minimum(a)
            n = sum.([diff(a; dims) for dims in 1:N])
            Z = norm(n)
            if Z != 0
                n /= Z
            end
        end
        P = n * n'
        inva = P * mean(1 ./ a) + (I - P) / mean(a)
        # reduce(vcat, [[inva[i, j] for j = 1:i] for i = 1:N])
    end
    [getindex.(inva, i, j) for i = 1:N, j = 1:N]
end