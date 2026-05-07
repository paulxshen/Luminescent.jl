const ϵPEC = 1337f0
const TOL = 1.0f-3
const ϵ0 = 8.854187817f-12
const PIXELS = 10
isPEC(x::Number) = abs(x) > ϵPEC / 5
isPEC(x) = isPEC.(x)
# const 𝐟 = BFloat16
const 𝐟 = Float32
const FF = Float32

const BREAK = "="^40

UnPack.unpack(x::AbstractDict, k::Val{f}) where f = x[f]