"""
    function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TMz (Ez, Hx, Hy) or :TEz (Hz, Ex, Ey)
"""
function setup(boundaries, sources, monitors, dx, sz0, polarization=nothing; ϵ=1, μ=1, σ=0, σm=0, F=Float32, Courant=0.5f0, kw...)# Courant number)
    # sz0 = round.(Int, L ./ dx)
    # sz0=size(ϵ)
    a = ones(F, sz0)
    # if isa(μ,Number+)
    μ *= a
    σ *= a
    σm *= a
    ϵ *= a
    d = length(sz0)
    if d == 1
        fk = (:Ez, :Hy)
    elseif d == 2
        if polarization == :TMz
            fk = (:Ez, :Hx, :Hy)
        elseif polarization == :TEz
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
    starts = Dict([k => ones(Int, d) for k = fk])
    sizes = Dict([k => collect(sz0) for k = fk])

    for b = boundaries
        for i = b.dims
            if typeof(b) == Periodic
                db[i, :] = [b, b]
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
                push!(geometry_padding[:σ], OutPad(ReplicateRamp(F(b.σ)), l, r,))
                push!(geometry_padding[:σm], OutPad(ReplicateRamp(F(b.σ)), l, r,))
                if j == 1
                    for k = keys(starts)
                        starts[k][i] += n
                    end
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
    for i = 1:d
        if d == 3
            l = OutPad(:replicate, Int.(1:d .== i), fill(0, d),)
            r = OutPad(:replicate, fill(0, d), Int.(1:d .== i),)
            p = [l, r]
            append!(geometry_padding[:μ], p)
            append!(geometry_padding[:σm], p)
            append!(geometry_padding[:σ], p)
            append!(geometry_padding[:ϵ], p)
            # push!(geometry_padding[:ϵ], p)
        end
    end
    # for (k, v) = pairs(geometry_padding)
    #     for p = v
    #         for v = values(starts)
    #             v .+= p.l .+ p.r
    #         end
    #     end
    # end
    for (k, v) = pairs(field_padding)
        for p = v
            sizes[k] += p.l .+ p.r
            starts[k] += p.l
        end
    end
    # Ie = [a:b for (a, b) = zip(start, start .+ esz .- 1)] # end
    # Ih = [a:b for (a, b) = zip(start, start .+ hsz .- 1)] # end



    starts[:Jx] = starts[:Ex]
    starts[:Jy] = starts[:Ey]
    starts[:Jz] = starts[:Ez]
    sizes[:Jx] = sizes[:Ex]
    sizes[:Jy] = sizes[:Ey]
    sizes[:Jz] = sizes[:Ez]
    sizes = NamedTuple([k => Tuple(sizes[k]) for (k) = keys(sizes)])
    fields = NamedTuple([k => zeros(F, Tuple(sizes[k])) for (k) = fk])
    u = collect(values(fields))
    if d == 1
        pf = u -> [u[1] .* u[2]]
    elseif d == 3
        u0 = [u[1:3], u[4:6]]
    else
        if polarization == :TMz
            u -> [-u[1] .* u[3], u[2] .* u[1]]
        else
            pf = u -> [u[1] .* u[3], -u[2] .* u[1]]
        end
    end
    geometry_splits = Dict{Symbol,Vector{Vector{UnitRange{Int64}}}}()
    function gs(sz, k)
        l = 1 .- sum(field_padding[k]) do p
            p.l
        end
        r = 1 .- sum(field_padding[k]) do p
            p.r
        end
        [ax[1+l:end-r] for (ax, l, r) = zip(Base.oneto.(sz), l, r)]
    end
    geometry_sizes = NamedTuple([k => sz0 .+ sum(geometry_padding[k]) do p
        p.l + p.r
    end for k = keys(geometry_padding)])
    geometry_splits[:μ] = geometry_splits[:σm] = gs.((geometry_sizes[:μ],), [:Hx, :Hy, :Hz])
    geometry_splits[:ϵ] = geometry_splits[:σ] = gs.((geometry_sizes[:ϵ],), [:Ex, :Ey, :Ez])
    source_effects = DefaultDict(() -> SourceEffect[])
    for s = sources
        for (k, v) = pairs(SourceEffect(s, dx, sizes, starts, sz0))
            append!(source_effects[k], v)
        end
    end
    c = 0
    save_info = Int[]

    monitor_instances = [
        begin

            idxs = Dict([
                k => map(starts[k], m.span) do s, x
                    x = round.(Int, x / dx)
                    if isa(x, Real)
                        s .+ x
                    else
                        s+x[1]:s+x[2]
                    end
                end for k = fk
            ])
            # idxs = NamedTuple([k => map(starts[k], m.span) do s, x
            #     x = round.(Int, x / dx)
            #     if isa(x, Real)
            #         s .+ x
            #     else
            #         s+x[1]:s+x[2]
            #     end
            # end for k = fk])
            centers = [k => round.(Int, mean.(v)) for (k, v) = pairs(idxs)]
            MonitorInstance(idxs, centers, m.normal, dx, m.label)
        end for m = monitors
    ]
    dt = dx * Courant
    sz0 = Tuple(sz0)
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
    (; μ, σ, σm, ϵ, geometry_padding, field_padding, geometry_splits, source_effects, monitor_instances, u0, save_info, fields, dx, dt, kw...)
end
function apply!(p; kw...)
    [apply!(p[k], kw[k]) for k = keys(kw)]
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function apply(p; kw...)
    [apply(p[k], kw[k]) for k = keys(kw)]
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function mark(p; kw...)
    [mark(p[k], kw[k]) for k = keys(kw)]
end
function mark(v, a)
    l = sum(v) do p
        p.l
    end
    r = sum(v) do p
        p.r
    end
    PaddedArray(a, l, r)
end

Base.strip(a::AbstractArray, l, r=l) = a[[ax[1+l:end-r] for (ax, l, r) = zip(axes(a), l, r)]...]
function apply(v, a)
    [a[i...] for i = v]
    # collect(getindex.((a,), zip(v...)...))
end