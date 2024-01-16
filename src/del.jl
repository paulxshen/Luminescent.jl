
VF = Union{AbstractVector{<:AbstractArray},Tuple{<:AbstractArray},NamedTuple}
struct Del
    Δ::AbstractArray
    cd::Bool

    function Del(Δ::AbstractArray, cd::Bool=false)
        new(Δ, cd)
    end
end
@functor Del

function pdiff(a::AbstractArray, ; dims, cd=false)
    if cd
    else
        diff(a; dims)
    end
end
function pdiff(a::PaddedArray, ; dims, cd=false)
    @unpack a, l, r = a
    I = [i == dims ? (:) : a:b for (i, (a, b)) = enumerate(zip(l .+ 1, size(a) .- r))]
    pdiff(a[I...]; dims, cd)
end
function (m::Del)(a::AbstractArray{<:Number}, p=*)
    n = length(m.Δ)
    if n == 1
        return pdiff(a) / m.Δ[1]
    end

    I = [ax[begin:end-1] for ax = axes(a)[1:n]]
    if p == *
        return [pdiff(a, dims=i) for i = 1:n] ./ m.Δ
    elseif p == cross
        dx, dy = m.Δ
    end
end

function (m::Del)(a, p=*)
    # function (m::Del)(a::AbstractVector{<:AbstractArray}, p=*)
    @unpack cd = m
    n = length(m.Δ)
    # I = [ax[begin:end-1] for ax = axes(first(a))]
    if p == dot
        return sum([pdiff(a[i], dims=i) for i = 1:n] ./ m.Δ)
    elseif p == cross
        if n == 2
            dx, dy = m.Δ
            if length(a) == 1
                return [pdiff(a[1], dims=2) / dy, -pdiff(a[1], dims=1) / dx]
            else
                u, v = a
                return [pdiff(v, dims=1) / dx - pdiff(u, dims=2) / dy]
                # return [pdiff(v, dims=1) / dx - pdiff(u, dims=2) / dy]
            end
        elseif n == 3
            dx, dy, dz = m.Δ
            u, v, w = a
            uy = pdiff(u, dims=2) / dy
            uz = pdiff(u, dims=3) / dz
            vx = pdiff(v, dims=1) / dx
            vz = pdiff(v, dims=3) / dz
            wx = pdiff(w, dims=1) / dx
            wy = pdiff(w, dims=2) / dy
            return [wy - vz, uz - wx, vx - uy]
        end
    end
end
function LinearAlgebra.cross(a::VF, b::VF)
    u, v, w = b
    x, y, z = a
    return [w .* y - v .* z, u .* z - w .* x, v .* x - u .* y]
end
function LinearAlgebra.dot(m::Del, a)
    m(a, dot)
end

function LinearAlgebra.cross(m::Del, a)
    m(a, cross)
end

"""
    Del(resolutions::AbstractVector)
    Del(cell::AbstractMatrix)

constructs ∇ Delerator (derivative, gradient, divergence, curl) using central pdifference stencil. Because the stencil is of length 3 in each dimension, the result is 2 elements shorter in each dimension than the input. To instead retain the same size, use `border=:smooth` which pads the input

# Example
## 1d derivative
```
dx = 0.1
x = 0:dx:.5
y = x .^ 2
d = Del(dx)
@test d(y)≈[ 0.2, 0.4,0.6,0.8]
@test d(y, border=:smooth) ≈ [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
```

## 2d gradient
```
dy = dx = 0.1
a = [x^2 + y^2 for x in 0:dx:0.5, y in 0:dy:0.5]
∇ = Del([dx, dy])
∇(a)

#=
4×4 Matrix{SVector{2, Float64}}:
 [0.2, 0.2]  [0.2, 0.4]  [0.2, 0.6]  [0.2, 0.8]
 [0.4, 0.2]  [0.4, 0.4]  [0.4, 0.6]  [0.4, 0.8]
 [0.6, 0.2]  [0.6, 0.4]  [0.6, 0.6]  [0.6, 0.8]
 [0.8, 0.2]  [0.8, 0.4]  [0.8, 0.6]  [0.8, 0.8]
=#
```

## 2d divergence, curl
curl of 2d Vector field is taken to be vorticity scalar field. In 3d curl produces vorticity vector field.
```
a = collect.([
    (0, 0) (-1, 0) (0, 0)
    (0, -1) (0, 0) (0, 1)
    (0, 0) (1, 0) (0, 0)
])
∇ = Del([1, 1])
@test ∇ ⋅ a ≈ ∇(a, dot) ≈ [2]'
@test ∇ × a ≈ ∇(a, cross) ≈ [0]'

a = collect.([
    (0, 0) (0, 1) (0, 0)
    (-1, 0) (0, 0) (1, 0)
    (0, 0) (0, -1) (0, 0)
])
∇ = Del([1, 1])
@test ∇ ⋅ a ≈ [0]' 
@test ∇ × a ≈ [-2]' 
```
"""

Base.:*(a::VF, b::VF) =
    broadcast(values(a), values(b)) do a, b
        a .* b
    end
Base.:*(a::VF, b::AbstractArray{<:Number}) = a * [b]
Base.:*(a::AbstractArray{<:Number}, b::VF) = [a] * b
Base.:/(a::VF, b::VF) =
    broadcast(values(a), values(b)) do a, b
        a ./ b
    end
Base.:/(a::VF, b::AbstractArray{<:Number}) = a / [b]
# Base./(a::AbstractArray{<:Number},b::VF) =[a]*b
*