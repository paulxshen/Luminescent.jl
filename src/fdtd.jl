"""
    function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TMz (Ez, Hx, Hy) or :TEz (Hz, Ex, Ey)
"""
function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32, Courant=0.5f0, kw...)# Courant number)
    sz0 = round.(Int, L ./ dx)
    d = length(sz0)
    esz = collect(sz0)
    hsz = copy(esz)
    npad = zeros(Int, d, 2)
    nodes = fill(:U, d, 2)
    db = Any[PML(j * i,) for i = 1:d, j = (-1, 1)]

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
                npad[i, j] += round(Int, b.d / dx)
                # hpad[i, j] .+= b.d
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

            hsz[i] += 1
            μ = pad(μ, :replicate, fill(0, d), Int.(1:d .== i))
        elseif nodes[i, :] == [:H, :H]

            esz[i] += 1
            ϵ = pad(ϵ, :replicate, fill(0, d), Int.(1:d .== i))
            σ = pad(σ, :replicate, fill(0, d), Int.(1:d .== i))
        end
    end

    start = npad[:, 1] .+ 1
    Ie = [a:b for (a, b) = zip(start, start .+ esz .- 1)] # end
    Ih = [a:b for (a, b) = zip(start, start .+ hsz .- 1)] # end

    esz0 = Tuple(esz)
    hsz0 = Tuple(hsz)
    esz += sum(npad, dims=2)
    hsz += sum(npad, dims=2)
    esz = Tuple(esz)
    hsz = Tuple(hsz)
    if d == 1
        fields = (;
            Ez=zeros(F, esz),
            Hy=zeros(F, hsz),
            # Jz=zeros(F, esz),
        )
    elseif d == 2
        if polarization == :TMz
            fields = (;
                Ez=zeros(F, esz),
                Hx=zeros(F, hsz),
                Hy=zeros(F, hsz),
                # Jz=zeros(F, esz),
            )
        elseif polarization == :TEz
            fields = (;
                Ex=zeros(F, esz),
                Ey=zeros(F, esz),
                Hz=zeros(F, hsz),
                # Jx=zeros(F, esz),
                # Jy=zeros(F, esz),
            )
        end
    else
        fields = (;
            Ex=zeros(F, esz),
            Ey=zeros(F, esz),
            Ez=zeros(F, esz),
            Hx=zeros(F, hsz), Hy=zeros(F, hsz), Hz=zeros(F, hsz),
            # Jx=zeros(F, esz),
            # Jy=zeros(F, esz),
            # Jz=zeros(F, esz),
        )
    end
    fields = NamedTuple([k => collect(v) for (k, v) = pairs(fields)])

    field_padding = DefaultDict(() -> Padding[])
    geometry_padding = DefaultDict(() -> Padding[])
    for i = 1:d
        xyz = para = perp = [:x, :y, :z]

        perp = [popat!(para, i)]
        for j = 1:2
            b = db[i, j]
            if isa(b, PML)
                n = round(Int, b.d / dx)
                l = j == 1 ? [i == a ? n : 0 for a = 1:d] : zeros(Int, d)
                r = j == 2 ? [i == a ? n : 0 for a = 1:d] : zeros(Int, d)
                info = lazy = false
                push!(geometry_padding[:ϵ], Padding(:ϵ, :replicate, l, r, info, lazy))
                push!(geometry_padding[:μ], Padding(:μ, :replicate, l, r, info, lazy))
                push!(geometry_padding[:σ], Padding(:σ, ReplicateRamp(F(b.σ)), l, r, info, lazy))
                push!(geometry_padding[:σm], Padding(:σm, ReplicateRamp(F(b.σ)), l, r, info, lazy))
            end
            l = j == 1 ? Int.((1:d) .== i) : zeros(Int, d)
            r = j == 2 ? Int.((1:d) .== i) : zeros(Int, d)
            t = typeof(b)


            f = nodes[i, j]
            for k = keys(fields)
                if startswith(String(k), String(f)) && k[2] in para
                    info = true
                    lazy = false
                    if t == Periodic
                        push!(field_padding[k], Padding(k, :periodic, l, r, info, lazy))
                    else
                        push!(field_padding[k], Padding(k, 0, l, r, info, lazy))
                    end
                end
            end
        end
    end
    stop = start .+ sz0 .- 1
    source_effects = [SourceEffect(s, dx, esz, start, stop) for s = sources]

    c = 0
    save_info = Int[]
    monitor_info = [
        map(start, m.span) do s, x
            idxs = s .+ round.(Int, x / dx)
            center = round.(Int, mean.(idxs))
            (; idxs, center)
        end for m = monitors
    ]
    dt = dx * Courant
    sz0 = Tuple(sz0)
    (; geometry_padding, field_padding, source_effects, monitor_info, save_info, fields, Ie, Ih, dx, esz0, hsz0, esz, hsz, dt, kw...)
end
