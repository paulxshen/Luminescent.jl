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


# power(m, u) = sum(sum(pf([u[m.idxs[k]...] for (u, k) = zip(u, fk)]) .* m.normal))
"""
    function power(m::MonitorInstance, u)

total power (Poynting flux) passing thru monitor surface
"""
function power(m::MonitorInstance, u)
    E = [u[1][i][m.idxs[k]...] for (i, k) = enumerate([:Ex, :Ey, :Ez])]
    H = [u[2][i][m.idxs[k]...] for (i, k) = enumerate([:Hx, :Hy, :Hz])]
    # E, H = u
    # E = [E[i][m.idxs[k]...] for (i, k) = enumerate([:Ex, :Ey, :Ez])]
    # H = [H[i][m.idxs[k]...] for (i, k) = enumerate([:Hx, :Hy, :Hz])]
    # E = [a[m.idxs[k]...] for (a, k) = zip(E, [:Ex, :Ey, :Ez])]
    # H = [a[m.idxs[k]...] for (a, k) = zip(H, [:Hx, :Hy, :Hz])]
    sum(sum(E × H .* m.normal))
end
# power(m, u) = sum(sum(pf([u[i...] for (u, i) = zip(u, collect(values(m.idxs)))]) .* m.normal))

function monitors_on_box(o, L)
    ox, oy, oz = o
    lx, ly, lz = L
    [
        Monitor([ox, [oy, oy + ly], [oz, oz + lz]], [-1, 0, 0]),
        Monitor([ox + lx, [oy, oy + ly], [oz, oz + lz]], [1, 0, 0]),
        Monitor([[ox, ox + lx], oy, [oz, oz + lz]], [0, -1, 0]),
        Monitor([[ox, ox + lx], oy + ly, [oz, oz + lz]], [0, 1, 0]),
        Monitor([[ox, ox + lx], [oy, oy + ly], oz], [0, 0, -1]),
        Monitor([[ox, ox + lx], [oy, oy + ly], oz + lz], [0, 0, 1]),
    ]
end