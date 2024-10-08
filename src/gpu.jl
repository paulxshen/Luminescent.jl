# gpu(F, x) = F(Flux.gpu(x))
# cpu(F, x) = F(Flux.cpu(x))
# using CUDA
# using GPUArraysCore
gpu(F, x) = cu(F.(x))
cpu(F, x) = Array(F.(x))
cpu(x::Number) = x
gpu(x::Number) = x
cpu(x::AbstractArray{<:Number}) = Array(x)
gpu(x::AbstractArray{<:Number}) = cu(x)
gpu(x::AbstractArray) = gpu.(x)
cpu(x::AbstractArray) = cpu.(x)


gpu(x) = isempty(propertynames(x)) ? x : fmap(gpu, x)
cpu(x) = isempty(propertynames(x)) ? x : fmap(cpu, x)
gpu(d::Map) = fmap(gpu, d)
cpu(d::Map) = fmap(cpu, d)
# gpu(v::AbstractVector{<:Collection}) = gpu.(v)
# cpu(v::AbstractVector{<:Collection}) = cpu.(v)

