T = AbstractArray{<:Number}
_cpu(x::T) = Array(x)
_gpu(x::T) = cu(x)
_gpu(x) = x
_cpu(x) = x
gpu(d) = fmap(_gpu, d, T)
cpu(d) = fmap(_cpu, d, T)

