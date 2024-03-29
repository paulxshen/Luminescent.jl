
function sandwich(mask, h1, h, h2, ϵ1, ϵ2)
    # a = ones(F, size(mask))
    a = 0 * mask .+ 1
    # hbox, hwg, hclad = h
    cat(repeat.(
            (a * ϵ1, mask, ϵ2 * a),
            1,
            1,
            (h1, h, h2)
        )..., dims=3)
end


_getindex(x::Real, a...) = x
_getindex(x, a...) = x[a...]
