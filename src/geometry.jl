struct SubpixelAveraging
    v
end

function apply(s::SubpixelAveraging, a::AbstractArray)
    d = ndims(a)
    # for (i, (do_smooth, pad_left, pad_right)) in enumerate(s.v)
    for (i, (l, r)) = enumerate(eachcol(s.v))
        select = i .== (1:d) |> Tuple
        if l == -1
            a = pad(a, :replicate, select, 0)
        end
        if r == 1
            a = pad(a, :replicate, 0, select)
        end
        if l != 0 && r != 0
            # a = conv(a, reshape([1, 1], 1 + select)) / 2
            a = (2selectdim(a, i, 1:(size(a, i)-1)) + diff(a, dims=i)) / 2
        end
    end
    a
end

function sandwich(layers, heights)
    sz = (1,)
    F = Float32
    for l = layers
        if isa(l, AbstractMatrix)
            sz = size(l)
            F = eltype(l)
        end
    end
    a = ones(F, sz)
    stack(vcat([repeat(isa(l, AbstractMatrix) ? l : l * a, n) for (l, n) = zip(layers, heights)]...))
end


# recursive_getindex(x::Real, a...) = x
# recursive_getindex(x, a...) = x[a...]
