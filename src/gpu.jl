Flux.gpu(d::Dictlike) = apply_func(Flux.gpu, d)
Flux.cpu(d::Dictlike) = apply_func(Flux.cpu, d)

S = Union{MonitorInstance,Collection}
Flux.gpu(v::AbstractVector{T}) where {T<:S} = gpu.(v)
Flux.cpu(v::AbstractVector{T}) where {T<:S} = cpu.(v)
function Flux.cpu(m::MonitorInstance)
    m.wavelength_modes = Flux.cpu(m.wavelength_modes)
    m
end
function Flux.gpu(m::MonitorInstance)
    m.wavelength_modes = Flux.gpu(m.wavelength_modes)
    m
end
