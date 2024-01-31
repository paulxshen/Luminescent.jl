using LinearAlgebra, UnPack
"""
    function Monitor(span, normal=nothing)

Constructs monitor which can span a point, line, surface, or volume

Args
- span
- normal: flux monitor direction
"""
struct Monitor
    span
    normal
    label
    function Monitor(span, normal=nothing, label="")
        new(span, normal, label)
    end
end
struct MonitorInstance
    idxs
    centers
    normal
    dx
    label
end