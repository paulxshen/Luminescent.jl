# gpu(F, x) = F(Flux.gpu(x))
# cpu(F, x) = F(Flux.cpu(x))
# using CUDA
# using GPUArraysCore
gpu(F, x) = cu(F.(x))
cpu(F, x) = Array(F.(x))

cpu(x::AbstractArray{<:Number}) = cpu(Float32, x)
gpu(x::AbstractArray{<:Number}) = gpu(Float32, x)
gpu(x::AbstractArray) = gpu.(x)
cpu(x::AbstractArray) = cpu.(x)


gpu(x) = isempty(propertynames(x)) ? x : fmap(gpu, x)
cpu(x) = isempty(propertynames(x)) ? x : fmap(cpu, x)
gpu(d::Dictlike) = fmap(gpu, d)
cpu(d::Dictlike) = fmap(cpu, d)
# gpu(v::AbstractVector{<:Collection}) = gpu.(v)
# cpu(v::AbstractVector{<:Collection}) = cpu.(v)

