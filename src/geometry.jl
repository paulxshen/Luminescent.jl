struct Staggering
    i
end

function apply(s::Staggering, a::AbstractArray)
    a[s.i...]
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
