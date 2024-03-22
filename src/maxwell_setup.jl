"""
    function maxwell_setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TM (Ez, Hx, Hy) or :TE (Hz, Ex, Ey)
"""
function maxwell_setup(boundaries, sources, monitors, dx, sz, polarization=:TE;
    ϵ=1, μ=1, σ=0, σm=0, F=Float32,
    ϵmin=1,
    # Courant=0.5,
    Courant=F(0.8√(ϵmin / length(sz))),# Courant number)
    kw...)
    dx = F(dx)

    a = ones(F, sz)
    # if isa(μ,Number+)
    μ *= a
    σ *= a
    σm *= a
    ϵ *= a
    d = length(sz)
    if d == 1
        fk = (:Ez, :Hy)
    elseif d == 2
        if polarization == :TM
            fk = (:Ez, :Hx, :Hy)
        elseif polarization == :TE
            fk = (:Ex, :Ey, :Hz)
        end
    else
        fk = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
    end
    Courant = F(Courant)

    global nodes = fill(:U, d, 2)
    db = Any[PML(j * i,) for i = 1:d, j = (-1, 1)]
    field_padding = DefaultDict(() -> InPad[])
    geometry_padding = DefaultDict(() -> OutPad[])
    fl = Dict([k => zeros(Int, d) for k = fk])
    flb = Dict([k => zeros(Int, d) for k = fk])
    frb = Dict([k => zeros(Int, d) for k = fk])
    lc = zeros(Int, d)
    sizes = Dict([k => collect(sz) for k = fk])

    for b = boundaries
        for i = b.dims
            if typeof(b) == Periodic
                db[i, :] = [Periodic(-abs(i)), Periodic(abs(i))]
            else
                if i > 0
                    db[i, 2] = b
                else
                    db[abs(i), 1] = b

                end
            end
        end
    end

    for i = 1:d
        for j = 1:2
            b = db[i, j]
            t = typeof(b)
            if t == PML
            elseif t == PEC
                nodes[i, j] = :E
            elseif t == PMC
                nodes[i, j] = :H
            end
        end
    end

    for i = 1:d
        if nodes[i, :] == [:U, :E]

            nodes[i, 1] = :H
        elseif nodes[i, :] == [:U, :H]

            nodes[i, 1] = :E
        elseif nodes[i, :] == [:E, :U]

            nodes[i, 2] = :H
        elseif nodes[i, :] == [:H, :U]

            nodes[i, 2] = :E
        elseif nodes[i, :] == [:U, :U]

            nodes[i, :] = [:E, :H]
        elseif nodes[i, :] == [:U, :U]

            nodes[i, :] = [:E, :H]
        elseif nodes[i, :] == [:E, :E]

        elseif nodes[i, :] == [:H, :H]

        end

    end



    for i = 1:d
        select = i .== 1:d
        xyz = para = perp = [:x, :y, :z]

        perp = [popat!(para, i)]
        for j = 1:2
            b = db[i, j]
            if isa(b, PML)
                n = round(Int, b.d / dx)
                l = j == 1 ? [i == a ? n : 0 for a = 1:d] : zeros(Int, d)
                r = j == 2 ? [i == a ? n : 0 for a = 1:d] : zeros(Int, d)
                info = lazy = false
                push!(geometry_padding[:ϵ], OutPad(:replicate, l, r,))
                push!(geometry_padding[:μ], OutPad(:replicate, l, r,))

                l1 = l .÷ 2
                r1 = r .÷ 2
                push!(geometry_padding[:σ], OutPad(ReplicateRamp(F(b.σ)), l1, r1,))
                push!(geometry_padding[:σm], OutPad(ReplicateRamp(F(b.σ)), l1, r1,))
                l2 = l - l1
                r2 = r - r1
                push!(geometry_padding[:σ], OutPad(F(b.σ), l2, r2,))
                push!(geometry_padding[:σm], OutPad(F(b.σ), l2, r2,))
                if j == 1
                    for k = keys(fl)
                        fl[k][i] += n
                    end
                    lc[i] += n
                end
                for k = keys(sizes)
                    sizes[k][i] += n
                end
            end
            l = j == 1 ? Int.((1:d) .== i) : zeros(Int, d)
            r = j == 2 ? Int.((1:d) .== i) : zeros(Int, d)
            t = typeof(b)


            f = nodes[i, j]
            for k = fk
                q = startswith(String(k), String(f))
                if (q ? k[2] in para : k[2] in perp)
                    if t == Periodic
                        lp = j == 1 ? select : zeros(Int, d)
                        rp = j == 2 ? select : zeros(Int, d)
                        push!(field_padding[k], InPad(:periodic, lp, rp,))
                    else
                        push!(field_padding[k], InPad(0, l, r,))
                    end
                end
            end
        end
    end

    p = OutPad(:replicate, ones(Int, d), ones(Int, d))
    push!(geometry_padding[:μ], p)
    push!(geometry_padding[:σm], p)
    push!(geometry_padding[:σ], p)
    push!(geometry_padding[:ϵ], p)

    for (k, v) = pairs(field_padding)
        for p = v
            sizes[k] += p.l .+ p.r
            fl[k] += p.l
            flb[k] += p.l
            frb[k] += p.r
        end
    end

    if d == 3
        fl[:Jx] = fl[:Ex]
        fl[:Jy] = fl[:Ey]
        fl[:Jz] = fl[:Ez]
        sizes[:Jx] = sizes[:Ex]
        sizes[:Jy] = sizes[:Ey]
        sizes[:Jz] = sizes[:Ez]
    elseif d == 2
        fl[:Jx] = fl[:Ex]
        fl[:Jy] = fl[:Ey]
        sizes[:Jx] = sizes[:Ex]
        sizes[:Jy] = sizes[:Ey]
    end

    sizes = NamedTuple([k => Tuple(sizes[k]) for (k) = keys(sizes)])
    fields = NamedTuple([k => zeros(F, Tuple(sizes[k])) for (k) = fk])
    if d == 1
        pf = u -> [u[1] .* u[2]]
    elseif d == 3
        u0 = dict([
            :E => dict([k => zeros(F, Tuple(sizes[k])) for k = (:Ex, :Ey, :Ez)]),
            :H => dict([k => zeros(F, Tuple(sizes[k])) for k = (:Hx, :Hy, :Hz)]),
            :J => dict([k => zeros(F, Tuple(sizes[k])) for k = (:Jx, :Jy, :Jz)]),
        ])
        # u0 = [u[1:3], u[4:6]]
    else
        if polarization == :TM
            u -> [-u[1] .* u[3], u[2] .* u[1]]
        else
            u0 = dict([
                :E => dict([k => zeros(F, Tuple(sizes[k])) for k = (:Ex, :Ey)]),
                :H => dict([k => zeros(F, Tuple(sizes[k])) for k = (:Hz,)]),
                :J => dict([k => zeros(F, Tuple(sizes[k])) for k = (:Jx, :Jy)]),
            ])
        end
    end

    geometry_staggering = Dict{Symbol,Any}()
    geometry_sizes = NamedTuple([k => sz .+ sum(geometry_padding[k]) do p
        p.l + p.r
    end for k = keys(geometry_padding)])

    geometry_staggering[:μ] =
        dict([Symbol("μ$k") => [ax[2-l:end-1+r] for (ax, l, r) = zip(Base.oneto.(geometry_sizes[:μ]), flb[k], frb[k])] for k = keys(u0[:H])])
    geometry_staggering[:σm] =
        dict([Symbol("σm$k") => [ax[2-l:end-1+r] for (ax, l, r) = zip(Base.oneto.(geometry_sizes[:μ]), flb[k], frb[k])] for k = keys(u0[:H])])
    geometry_staggering[:ϵ] =
        dict([Symbol("ϵ$k") => [ax[2-l:end-1+r] for (ax, l, r) = zip(Base.oneto.(geometry_sizes[:ϵ]), flb[k], frb[k])] for k = keys(u0[:E])])
    geometry_staggering[:σ] =
        dict([Symbol("σ$k") => [ax[2-l:end-1+r] for (ax, l, r) = zip(Base.oneto.(geometry_sizes[:ϵ]), flb[k], frb[k])] for k = keys(u0[:E])])

    source_instances = SourceInstance.(sources, dx, (sizes,), (lc,), (fl,), (sz,); F)
    monitor_instances = MonitorInstance.(monitors, dx, (lc,), (flb,), (fl,); F)
    roi = MonitorInstance(Monitor(zeros(d), zeros(d), dx * sz), dx, lc, flb, fl; F)

    dt = dx * Courant
    sz = Tuple(sz)
    for (k, v) = pairs(field_padding)
        a = ones(F, size(fields[k]))
        i = filter(eachindex(v)) do i
            v[i].b == 0
        end
        i_ = sort(setdiff(eachindex(v), i))
        if !isempty(i)

            for i = i
                @unpack b, l, r = v[i]
                pad!(a, b, l, r)
            end
            field_padding[k] = [
                InPad(
                    0,
                    sum(getproperty.(v[i], :l)),
                    sum(getproperty.(v[i], :r)),
                    a),
                v[i_]...
            ]
        end
    end
    geometry_padding = NamedTuple(geometry_padding)
    field_padding = NamedTuple(field_padding)
    #    verbose &&
    @info """
 ====
 FDTD configs
 
 Lengths in characterstic wavelengths, times in characterstic periods, unless otherwise specified

 dx: $(dx|>d2) 
 dt: $(dt|>d2) 
 Courant number: $(Courant|>d2)

 Original array size of all fields in pixels: $sz
 Padded field array sizes in pixels:
 $sizes
 
 Boundaries:
 $(join("- ".*string.(db),"\n"))
 
 Sources:
 $(join("- ".*string.(sources),"\n"))
 
 Monitors:
 $(join("- ".*string.(monitors),"\n"))
 
 $footer
 ====
 """

    (; μ, σ, σm, ϵ,
        geometry_padding, field_padding, geometry_staggering,
        source_instances, monitor_instances,
        roi, u0, fields, dx, dt, kw...)
end

