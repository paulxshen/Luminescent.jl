using LinearAlgebra, UnPack
"""
    function Monitor(c, L, n=nothing)

Constructs monitor which can span a point, line, surface, or volume monitoring fields or power

Args
- c: origin or center of monitor
- L: physical dimensions of monitor
- n: flux monitor direction (eg normal to flux surface)
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
    function Monitor(c, L, n=nothing, label="")
        new(c, -L / 2, L / 2, n, label)
    end
end
struct MonitorInstance
    idxs
    centers
    n
    dx
    A
    label
end


# power(m, u) = sum(sum(pf([u[m.idxs[k]...] for (u, k) = zip(u, fk)]) .* m.n))
"""
    function power_density(m::MonitorInstance, u)

 power density (avg Poynting flux) passing thru monitor surface
"""
function power_density(m::MonitorInstance, u)
    E = [u[1][i][m.idxs[k]...] for (i, k) = enumerate([:Ex, :Ey, :Ez])]
    H = [u[2][i][m.idxs[k]...] for (i, k) = enumerate([:Hx, :Hy, :Hz])]
    # E, H = u
    # E = [E[i][m.idxs[k]...] for (i, k) = enumerate([:Ex, :Ey, :Ez])]
    # H = [H[i][m.idxs[k]...] for (i, k) = enumerate([:Hx, :Hy, :Hz])]
    # E = [a[m.idxs[k]...] for (a, k) = zip(E, [:Ex, :Ey, :Ez])]
    # H = [a[m.idxs[k]...] for (a, k) = zip(H, [:Hx, :Hy, :Hz])]
    mean(sum(E × H .* m.n))
end
"""
    function power(m::MonitorInstance, u)

total power (Poynting flux) passing thru monitor surface
"""
function power(m::MonitorInstance, u)
    m.A * power_density(m, u)
end
# power(m, u) = sum(sum(pf([u[i...] for (u, i) = zip(u, collect(values(m.idxs)))]) .* m.n))

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