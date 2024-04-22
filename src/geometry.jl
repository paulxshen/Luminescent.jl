
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


_getindex(x::Real, a...) = x
_getindex(x, a...) = x[a...]
