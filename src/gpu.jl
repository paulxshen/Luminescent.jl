
unsafe_free!(a) = 0
unsafe_free!(a::CUDA.CuArray) = CUDA.unsafe_free!(a)