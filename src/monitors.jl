using LinearAlgebra, UnPack
# For example, `a::Array[m::MonitorInstance, :Ex]` queries Ex on monitor m.
"""
    function Monitor(c, L; normal=nothing, label="")
    function Monitor(c, lb, ub; normal=nothing, label="")

Constructs monitor which can span a point, line, surface, or volume monitoring fields or power. 

Args
- c: origin or center of monitor
- L: physical dimensions of monitor
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c

- normal: flux monitor direction (eg normal to flux surface)
"""
struct Monitor
    c
    lb
    ub
    n
    label
    # function Monitor(c,, n=nothing, label="")
    #     new(c,-L/2,L/2, n, label)
    # end
    function Monitor(c, L; normal=nothing, label="")
        new(c, -L / 2, L / 2, normal, label)
    end
    function Monitor(c, lb, ub; normal=nothing, label="")
        new(c, lb, ub, normal, label)
    end
end
struct MonitorInstance
    i
    c
    fi
    n
    dx
    A
    label
end
@functor MonitorInstance (n,)
Base.getindex(a::Union{AbstractArray,GPUArraysCore.AbstractGPUArray}, m::MonitorInstance) = a[m.i...]
Base.getindex(a::GPUArraysCore.AbstractGPUArray, m::MonitorInstance) = a[m.i...]
Base.getindex(a, m::MonitorInstance, k) = a[m.fi[k]...]

# power(m, u) = sum(sum(pf([u[m.i[k]...] for (u, k) = zip(u, fk)]) .* m.n))
"""
    function power_density(m::MonitorInstance, u)

 power density (avg Poynting flux) passing thru monitor surface
"""
function power_density(m::MonitorInstance, u)
    E = [u[1][i][m.fi[k]...] for (i, k) = enumerate([:Ex, :Ey, :Ez])]
    H = [u[2][i][m.fi[k]...] for (i, k) = enumerate([:Hx, :Hy, :Hz])]

    # E = getindex.(u[1], Ref.(m.i)...)
    # H = getindex.(u[2], Ref.(m.i)...)
    # E, H = u

    r = mean(sum((E × H) .* m.n))
    r
end
"""
    function power(m::MonitorInstance, u)

total power (Poynting flux) passing thru monitor surface
"""
function power(m::MonitorInstance, u)
    m.A * power_density(m, u)
end
# power(m, u) = sum(sum(pf([u[i...] for (u, i) = zip(u, collect(values(m.i)))]) .* m.n))

function monitors_on_box(c, L)
    ox, oy, oz = c
    lx, ly, lz = L
    rx, ry, rz = L / 2
    [
        Monitor([ox - rx, oy, oz], [0, ly, lz], [-1, 0, 0]),
        Monitor([ox + rx, oy, oz], [0, ly, lz], [1, 0, 0]),
        Monitor([ox, oy - ry, oz], [lx, 0, lz], [0, -1, 0]),
        Monitor([ox, oy + ry, oz], [lx, 0, lz], [0, 1, 0]),
        Monitor([ox, oy, oz - rz], [lx, ly, 0,], [0, 0, -1,]),
        Monitor([ox, oy, oz + rz], [lx, ly, 0,], [0, 0, 1,]),
    ]
end