for T = (:Map, :(AbstractVector{<:AbstractArray}))
    @eval Base.:*(A::AbstractMatrix{<:AbstractArray}, v::$T) = [sum([a .* v for (a, v) = zip(r, values(v))]) for r in eachrow(A)]
end
Base.:*(A::AbstractMatrix{<:AbstractArray}, B::AbstractMatrix{<:AbstractArray}) = [sum([a .* b for (a, b) = zip(r, c)]) for r = eachrow(A), c = eachcol(B)]

LinearAlgebra.norm(a::AbstractArray{<:AbstractArray}) = sqrt.(sum(a) do a
    a .⋅ a
end)
