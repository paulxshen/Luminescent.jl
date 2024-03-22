
function sandwich(base, h, ϵ)
    # a = ones(F, size(base))
    a = 0 * base .+ 1
    ϵbox, ϵcore, ϵclad, = ϵ
    # hbox, hwg, hclad = h
    cat(repeat.(
            (a * ϵbox, base * ϵcore + (1 .- base) * ϵclad, ϵclad * a),
            1,
            1,
            h
        )..., dims=3)
end


_getindex(x::Real, a...) = x
_getindex(x, a...) = x[a...]
