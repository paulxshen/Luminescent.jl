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
        if l == -1 == r
            # a = conv(a, reshape([1, 1], 1 + select)) / 2
            # a = (2selectdim(a, i, 1:(size(a, i)-1)) + diff(a, dims=i)) / 2
            a = (2selectdim(a, i, 2:(size(a, i))) - diff(a, dims=i)) / 2
        elseif l == 1 == r
            a = (2selectdim(a, i, 1:(size(a, i)-1)) + diff(a, dims=i)) / 2
        end
    end
    a
end

