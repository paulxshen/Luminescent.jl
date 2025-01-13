using Zygote, CUDA, ChainRulesCore, Porcupine
function ChainRulesCore.rrule(::typeof(CUDA.fill), v, sz)
    y = CUDA.fill(v, sz)
    function pb(ȳ)
        NoTangent(), sum(ȳ), NoTangent()
    end
    return y, pb
end
Zygote.gradient(2.0) do a
    # sum(a * CUDA.fill(a, (2, 2)))
    sum(a * CUDA.fill(a, (2, 2)))
end