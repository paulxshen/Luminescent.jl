"""
    function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TMz (Ez, Hx, Hy) or :TEz (Hz, Ex, Ey)
"""
function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)
    sz = round.(Int, L ./ dx)
    d = length(sz)
    esz = collect(sz)
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

    esz += sum(npad, dims=2)
    hsz += sum(npad, dims=2)
    esz = Tuple(esz)
    hsz = Tuple(hsz)
    if d == 1
        fields = (;
            Ez=zeros(F, esz),
            Hy=zeros(F, hsz),
            Jz=zeros(F, esz),)
    elseif d == 2
        if polarization == :TMz
            fields = (;
                Ez=zeros(F, esz),
                Hx=zeros(F, hsz),
                Hy=zeros(F, hsz),
                Jz=zeros(F, esz),)
        elseif polarization == :TEz
            fields = (;
                Ex=zeros(F, esz),
                Ey=zeros(F, esz),
                Hz=zeros(F, hsz),
                Jx=zeros(F, esz),
                Jy=zeros(F, esz),
            )
        end
    else
        fields = (;
            Ex=zeros(F, esz),
            Ey=zeros(F, esz),
            Ez=zeros(F, esz),
            Hx=zeros(F, hsz), Hy=zeros(F, hsz), Hz=zeros(F, hsz),
            Jx=zeros(F, esz),
            Jy=zeros(F, esz),
            Jz=zeros(F, esz),)
    end
    fields = NamedTuple([k => collect(v) for (k, v) = pairs(fields)])

    field_padding = Padding[]
    geometry_padding = Padding[]
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
                push!(geometry_padding, Padding(:ϵ, :replicate, l, r, info, lazy))
                push!(geometry_padding, Padding(:μ, :replicate, l, r, info, lazy))
                push!(geometry_padding, Padding(:σ, ReplicateRamp(F(b.σ)), l, r, info, lazy))
                push!(geometry_padding, Padding(:σm, ReplicateRamp(F(b.σ)), l, r, info, lazy))
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
                        push!(field_padding, Padding(k, :periodic, l, r, info, lazy))
                    else
                        push!(field_padding, Padding(k, 0, l, r, info, lazy))
                    end
                end
            end
        end
    end
    source_effects = [SourceEffect(s, dx, sz, start) for s = sources]

    c = 0
    save_idxs = Int[]
    monitor_configs = MonitorConfig[]

    for m = monitors
        @unpack span, k = m
        sz = length.(filter(x -> !isa(x, Number), span))
        if isempty(sz)
            sz = Int[]
        end

        idxs = map(start, span) do s, x
            x = round.(Int, x / dx)
            if length(x) == 1
                s + x
            end
        end

        # if isa(fields, NamedTuple)
        #     fi = 0
        # else

        # ki = NamedTuple([k =>
        #     begin
        #         fa = fields[k]
        #         # fcinds = label2index(fields, "$f")
        #         linds = LinearIndices(fa)
        #         fcinds[linds[(idxs)...]]
        #     end
        #                  for (k,start,stop)=k
        # ])
        # end
        ki = findfirst.(isequal.(k), (keys(fields),))
        push!(monitor_configs, MonitorConfig(idxs, ki))
    end
    (; geometry_padding, field_padding, source_effects, monitor_configs, save_idxs, fields, Ie, Ih)
end
