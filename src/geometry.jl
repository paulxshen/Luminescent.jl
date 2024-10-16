
# function _apply_subpixel_averaging(lr, a::AbstractArray)
#     N = ndims(a)
#     lr = cpu(lr)
#     l, r = eachcol(lr)
#     a = pad(a, :replicate, 1)
#     _l = l + 1.5
#     _r = _l + r - l
#     getindexf(a, range.(_l, _r + 0.001)...)
# end
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


function apply_subpixel_averaging(geometry, fieldlims)
    namedtuple([k => begin
        global adsf = fieldlims
        lrs = ignore_derivatives() do
            E = r"E.*"
            H = r"H.*"
            f = (; ϵ=E, σ=E, μ=H, m=H, invϵ=E)[k]
            fieldlims(f)
        end
        @assert length(lrs) > 0
        a = geometry[k]
        map(lrs) do lr
            crop(a, lr[:, 1] - 0.5, size(a) - 0.5 - lr[:, 2])
        end
        # end for (k, a) = pairs(geometry)])
    end for k = keys(geometry)])
end

function pad_geometry(geometry, geometry_padvals, geometry_padamts, ratio=1)
    namedtuple([
        k => begin
            if k in keys(geometry_padvals) && k in keys(geometry_padamts)
                pad(geometry[k], geometry_padvals[k][1], eachcol(geometry_padamts[k] * ratio)...)
            else
                geometry[k]
            end
        end for k = keys(geometry)])
end

function _apply_field_padvals(p::AbstractVector{<:InPad}, a::AbstractArray; nonzero_only=false)
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

function apply_field_padvals(fps, fs; kw...)
    # namedtuple([k => apply(fps[k], fs[k]; kw...) for k = keys(fs)])
    dict([k => k in keys(fps) ? _apply_field_padvals(fps[k], fs[k]; kw...) : fs[k]
          for k = keys(fs)])
end

function tensorinv(a, ratio, fieldlims)
    N = ndims(a)
    # T = typeof(a)
    F = eltype(a)
    # sz = size(a) .÷ ratio
    # margin = ratio ÷ 2
    # A = pad(a, :replicate, margin)
    # A = cpu(A)

    # v = [
    #     begin
    #         lr = cpu(lr)
    #         l, r = eachcol(lr)
    #         start = Int.((l - 0.5) * ratio) + margin + 1
    #         finish = size(A) - margin + Int.((r - sz + 0.5) * ratio)
    #         # @show start, finish, size(A), margin, l, r, sz
    #         _a = A[range.(start, finish)...]
    #         downsample(_a, ratio) do a
    #             n = ignore_derivatives() do
    #                 if maximum(a) == minimum(a)
    #                     zeros(F, N)
    #                 else
    #                     n = mean.([diff(a; dims) for dims in 1:N])
    #                     Z = norm(n)
    #                     # if Z != 0
    #                     if Z > 1e-2
    #                         n / Z
    #                     else
    #                         zeros(F, N)
    #                     end
    #                 end
    #             end
    #             P = n * n'
    #             P * mean(1 ./ a) + (I - P) / mean(a)
    #         end
    #     end for lr = fieldlims(r"E.*")
    # ]
    # T.([getindex.(v[j], i, j) for i = 1:N, j = 1:N])
    @assert !any(iszero, a)
    @assert !any(isnan, a)
    println(minimum(a))
    return [i == j ? 1 ./ a : zeros(F, size(a)) for i = 1:N, j = 1:N]
    # [i == j ? 1 ./ downsample(a, ratio) : zeros(F, size(a) .÷ ratio) for i = 1:N, j = 1:N]
end