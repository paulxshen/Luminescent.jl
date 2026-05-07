cummin(v, n) = minimum(mean.([v[i:i+n-1] for i = 1:length(v)-n+1]))
struct Deconv
    A
    ratio
    function Deconv(fs::AbstractArray{<:F}, tcut; ratio=2) where F
        _fs = resize(fs, size(fs) * ratio; inbound=true)
        _nf = length(_fs)

        dfs = pad(diff(fs), :replicate, 1)
        dfs = (dfs[1:end-1] + dfs[2:end]) / 2
        df = mean(dfs)
        nf = length(fs)

        f1 = 1 / (tcut)
        flim = 3f1 |> F

        nconv = round(Int, 2flim / df) + 1
        @debug (; nconv, nf, flim)
        # if 8nconv > nf
        #     println("T too small for time extrapolation which won't be performed")
        #     return LinearAlgebra.I
        # end
        n = nf * nconv

        I = zeros(Int, n)
        J = zeros(Int, n)
        V = zeros(complex(F), n)
        fmin, fmax = extrema(fs)

        c = 1
        for (i, _f) = enumerate(_fs)
            _flim = min(flim, _f - fmin, fmax - _f) + df / 1000
            @assert _flim > 0
            a = searchsortedfirst(fs, _f - _flim)
            b = searchsortedfirst(fs, _f + _flim) - 1
            @assert a <= b
            # b = max(a, b)

            v = fs[a:b] - _f
            h = dfs[a:b] .* sinc.(v / f1)
            h /= sum(h)
            h = cispi.(v * tcut) .* h
            for (j, w) = zip(a:b, h)
                I[c] = i
                J[c] = j
                V[c] = w
                c += 1
            end
        end

        I = I[1:(c-1)]
        J = J[1:(c-1)]
        V = V[1:(c-1)]
        I = _nf - I + 1
        J = nf - J + 1
        new(sparse(I, J, V), ratio)
    end
end
(m::Deconv)(v) = m.A \ resize(v, size(v) * m.ratio; inbound=true)

function setup(λ, λs, bbox, nres, geometry, boundaries, sources, monitors, modes=[], ;
    approx_2D_mode=nothing, z=nothing, name,
    tmax=nothing, energy_decay_threshold=nothing, path_length_multiple=nothing,
    F=Float32,
    relative_courant=0.9,
    array=Array,
    # subpixel_smoothing
    relative_pml_depths=1, path,
    verbosity=2, saveat=nothing, checkat=1000,
    bg=(ϵ=1, μ=1, σ=0, m=0), subpixel_smoothing=true,
    time_extrapolation=true,
    gradckptat=nothing,
    designs=[], targets=[],
    timestamp=nothing,
)
    ENV["start"] = string(time())
    # check_ping()

    bg = namedtuple(vmap(FF, bg))

    N = size(bbox, 1)

    L = diff(bbox, dims=2)
    L /= λ
    # FF = BFloat16
    println("""
    $BREAK
    setting up simulation...
    """)

    MODES = joinpath(path, "modes")

    boundaries = permutedims(reduce(hcat, map(boundaries) do x
            Symbol.(x isa Union{Tuple,AbstractArray} ? collect(x) : [x, x])
        end), (2, 1))

    meshes = Scale(1 / λ).(getindex.(geometry, 1))
    meshprops = getindex.(geometry, 2)

    npoles = maximum(meshprops, init=0) do d
        v = d(:dispersion)
        isnothing(v) ? 0 : length(v)
    end
    meshprops = map([meshprops..., bg]) do d
        d = kvmap(d) do k, v
            if v isa AbstractArray
                if v[1] isa AbstractArray
                    v = stack(v)
                end
            end
            Symbol(k) => FF(v)
        end
        if npoles > 0
            v = d[:dispersion]
            isnothing(v) && (v = [])
            v = vcat(collect.(v), fill([0.0f0, 0.0f0, nothing], npoles - length(v)))
            for v = v
                v[2] = v[2]^2
            end
            d[:dispersion] = v
        end
        d |> namedtuple
    end

    n = getindex.(meshprops, :mesh_density)
    nmax = maximum(maximum.(n))
    bbox /= λ
    z /= λ

    bbox, nres, relative_courant, relative_pml_depths = FF.((bbox, nres, relative_courant, relative_pml_depths))

    tol = 1 / nres / nmax
    for m = meshes
        x = boundingbox(m)
        # @debug (; x)
        a = ustrip.(getfield(coords(x.min), :coords))
        b = ustrip.(getfield(coords(x.max), :coords))
        tol = min(tol, minimum(b - a))
    end
    tol /= 100

    rulers = makemesh(zip(meshes, n), bbox |> FF, nres |> FF; bg=bg(:mesh_density, √(bg(:ϵ))), tol) .|> FF
    # @debug rulers * λ
    @debug length.(rulers)
    centers = centroids.(rulers)
    deltas = diff.(rulers)
    @debug extrema.(deltas) * λ

    ratio = 4
    _deltas = map(deltas) do v
        reduce(vcat, fill.(v / ratio, ratio))
    end
    _rulers = cumsum.(vcat.(first.(rulers), _deltas))

    # weights = centroids(pad(prod.(Base.product(deltas...)), :replicate, 1); dims=1:N)
    # @debug deltas[end] * λ
    sz = length.(rulers) - 1
    sz0 = Tuple(sz)
    if !isnothing(approx_2D_mode)
        approx_2D_mode = Symbol(approx_2D_mode)
    end

    if N == 1
        field_names = (:Ez, :Hy)
    elseif N == 2
        if approx_2D_mode == :TM
            Enames = (:Ez,)
            Hnames = (:Hx, :Hy)
            field_names = (:Ez, :Hx, :Hy)
        elseif approx_2D_mode == :TE
            Enames = (:Ex, :Ey)
            Hnames = (:Hz,)
            field_names = (:Ex, :Ey, :Hz)
        end
    else
        Enames = (:Ex, :Ey, :Ez)
        Hnames = (:Hx, :Hy, :Hz)
        field_names = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
    end

    nodes = fill(:U, N, 2)
    # db = Any[PML(j * i, lpml[i, int(j / 2 + 1.5)], σpml, mpml) for i = 1:N, j = (-1, 1)]

    nothingarray = () -> Array{Any,2}(fill(nothing, N, 2))
    zeroarray = () -> zeros(Int, N, 2)

    offsets = OrderedDict{Symbol,Any}()
    boundvals = OrderedDict{Symbol,Any}()
    geometry_names = [:ϵ, :μ, :σ, :m, :invϵ, :invμ, :χ2, :χ3]
    # npoles > 0 && push!(geometry_names, :Cn)

    #
    padvals = OrderedDict{Symbol,Any}([
        :default => nothingarray(),
        map(geometry_names) do k
            k => nothingarray()
        end...
    ])
    # padamts = OrderedDict{Symbol,Any}([
    #     :default => zeroarray(),
    #     map(geometry_names) do k
    #         k => zeroarray()
    #     end...
    # ])
    padamts = zeroarray()

    if npoles > 0
        padvals[:Cn] = nothingarray()
        # padamts[:Cn] = zeroarray()
    end

    if N == 1
    elseif N == 3
        all_field_names = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
    else
        if approx_2D_mode == :TM
            all_field_names = (:Ez, :Hx, :Hy)
        else
            all_field_names = (:Ex, :Ey, :Hz)
        end
    end

    walls = deepcopy(boundaries)
    for (i, (a, b)) = enumerate(eachrow(walls))
        if a == b == :PML
            walls[i, :] = [:PMC, :PEC]
        elseif a == :PML
            if b == :PEC
                walls[i, 1] = :PMC
            elseif b == :PMC
                walls[i, 1] = :PEC
            else
                error("Invalid boundary condition combination")
            end
        elseif b == :PML
            if a == :PEC
                walls[i, 2] = :PMC
            elseif a == :PMC
                walls[i, 2] = :PEC
            else
                error("Invalid boundary condition combination")
            end
        end
    end

    xyz = "xyz"
    for k = (:E, :H)
        for i = 1:N
            a = Symbol("$k$(xyz[i])")
            if a in all_field_names
                v = map(1:N) do j
                    if iseven((k == :E) + (i == j))
                        b = walls[j, 1]
                        if b == :PEC
                            b = :symmetric
                        elseif b == :PMC
                            b = :antisymmetric
                        else
                            @assert b == :periodic
                        end
                        0, [b, nothing]
                    else
                        b = walls[j, 2]
                        if b == :PEC
                            b = 0
                        elseif b == :PMC
                            # b = :taper
                            b = :replicate
                        else
                            @assert b == :periodic
                        end
                        -0.5, [nothing, b]
                    end
                end
                offsets[a] = first.(v) .|> FF
                boundvals[a] = permutedims(reduce(hcat, last.(v)), (2, 1))
            end
        end
    end
    add_current_keys!(offsets)


    geometry = OrderedDict{Symbol,Any}()
    centgeom = OrderedDict{Symbol,Any}()
    println("meshing geometry")

    I = findall(isnothing, getindex.(meshprops, :ϵ))
    deleteat!(meshprops, I)
    deleteat!(meshes, I)

    n = getindex.(meshprops, :mesh_density)
    REL_SAMPLE_DIST = 0.15
    searchers = Searcher.(meshes, (REL_SAMPLE_DIST ./ n[1:end-1] / nres) |> FF, tol)
    bboxes = getfield.(searchers, :bbox) |> stack
    ixs = UInt8.(eachindex(meshes))
    ixbg = UInt8(length(meshes) + 1)
    funnel = RangeFinder(bboxes)

    println("assigning materials")
    @time "enums" Enums = supersamplemesh(centers, deltas, funnel, searchers, bbox, [ixs..., ixbg], offsets(r"E(.)"), false) #nc
    @time "enums" centenum = supersamplemesh(centers, deltas, funnel, searchers, bbox, [ixs..., ixbg], [0], false)[1] #nc
    # @time "_enums" _enums = supersamplemesh(funnel, searchers, ixs, _rulers; z, bg=ixbg)[1] #nc

    # nadj = zeros(Int, N, 2)
    # for i = 1:N
    #     for j = 1:2
    #         lext = maximum(sources) do s
    #             max(0, j == 1 ? bbox[i, 1] - s.origin[i] / λ : s.origin[i] - bbox[i, 2] / λ)
    #         end
    #     end
    # end

    for k = (:ϵ, :μ, :σ, :χ2, :χ3)
        v = getindex.(meshprops, k)
        v = [isa(v, Number) ? v : mean(diag(v)) for v in v]
        vmin, vmax = extrema(v)
        @debug "$k: $vmin to $vmax"
        if vmin == vmax
            geometry[k] = centgeom[k] = vmin
        else
            geometry[k] = [getindex.((v,), a) for a = Enums]
            centgeom[k] = [getindex.((v,), centenum)]
        end
    end

    if npoles == 0
        geometry[:Cn] = 0
    else
        geometry[:Cn] = map(1:npoles) do i
            v = getindex.(getindex.(meshprops, :dispersion), i)
            v = [[getindex.((getindex.(v, j),), a) for a in Enums] for j = 1:3]

            v[3] = map(v[3]) do a
                a = filter(!isnothing, a)
                isempty(a) && return false
                m = mode(a)
                ifelse.(isnothing.(a), (m,), a)
            end

            map(v) do v
                s = map(v) do a
                    p, q = extrema(a)
                    p == q && return p
                    a
                end
                if all(isa.(s, Number))
                    s[1]
                else
                    v
                end
            end
        end
    end
    # global _c = geometry[:Cn]
    # error("debug Cn")

    geometry[:m] = [zeros(FF, sz0)]
    v = geometry[:σ]
    if v isa Number
        geometry[:σ] = [fill(v, sz0)]
    end
    geometry[:invμ] = 1 ⊘ geometry(:μ)

    ϵ = isoval.(getindex.(meshprops, :ϵ))
    hasPEC = any(isPEC, ϵ)
    ϵmaxdielectric = maximum(abs, ϵ .* .!(isPEC(ϵ)))

    if subpixel_smoothing
        @time geometry[:invϵ] = supersamplemesh(centers, deltas, funnel, searchers, bbox, ϵ, offsets(r"E(.)"), true)
    else
        geometry[:invϵ] = 1 ⊘ geometry(:ϵ)
    end


    ϵmin, ϵmax = extrema(ϵ)
    μmin, μmax = extrema(getindex.(meshprops, :μ))
    # nmax = sqrt(ϵmaxcapped * μmax)
    nmax_dielectric = sqrt(ϵmaxdielectric * μmax)
    nmin = sqrt(ϵmin * μmin)
    @debug (; ϵmin, ϵmax, μmin, μmax, nmin, nmax_dielectric)

    dt = FF(nmin) / √(sum(v -> 1 / cummin(v, 1)^2, FF(deltas))) * FF(relative_courant)
    @assert dt > 0
    dt = 1 / ceil(1 / dt) |> FF

    flush = fill(false, N, 2)
    for m = monitors
        for i = 1:N
            v = m.origin[i] / λ
            for j = 1:2
                d = j == 1 ? deltas[i][1] : deltas[i][end]
                lim = j == 1 ? rulers[i][1] : rulers[i][end]
                if abs(v - lim) < d
                    flush[i, j] = true
                end
            end
        end
    end
    @debug flush

    𝓁 = 0.25 |> FF
    lpmlmin = maximum(λs / λ) / nmin * 𝓁
    dpmlmax = minimum(λs / λ) / nmax_dielectric / nres |> FF
    npmlmin = 6

    k0 = -log(1.0f-2) / (2𝓁)
    k = min(0.1 / dt, k0) |> FF
    lpmlmin *= k0 / k
    σpml = ϵmin * k
    mpml = μmin * k
    @debug (; σpml, mpml, nmin, nmax_dielectric, ϵmin, μmin, dt, lpmlmin, dpmlmax, npmlmin)

    @debug (; z)
    for i = 1:N
        select = i .== 1:N
        xyz = para = perp = [:x, :y, :z]
        perp = [popat!(para, i)]
        for j = 1:2
            b = boundaries[i, j]
            relpml = relative_pml_depths(i)(j)
            if b == :PML
                nzero = 1
                d = j == 1 ? deltas[i][1] : deltas[i][end]
                v = [d]
                dmax = max(d, dpmlmax)
                while length(v) < npmlmin * relpml || sum(v) < (lpmlmin) * relpml
                    # if sum(v) - v[1] < lext * relpml
                    #     nzero += 1
                    #     r = 1
                    # else
                    #     r = 1.5
                    # end
                    push!(v, min(dmax, 1.5v[end]))
                end
                # nadj[i, j] = nzero - 1

                npml = length(v)
                v = cumsum(v)

                if j == 1
                    prepend!(rulers[i], rulers[i][1] - reverse(v))
                else
                    append!(rulers[i], rulers[i][end] + v)
                end

                padamts[i, j] = npml
                for k = geometry_names
                    # padamts[k][i, j] = npml
                    if k ∉ (:σ, :m)
                        padvals[k][i, j] = :replicate
                    end
                end

                if npoles > 0
                    # padamts[:Cn][i, j] = npml
                    padvals[:Cn][i, j] = :replicate
                end

                @debug (; nzero, npml)
                padvals[:σ][i, j] = Ramp(σpml * relpml, nzero)
                padvals[:m][i, j] = Ramp(mpml * relpml, nzero)

                sz[i] += npml
            end
        end
    end
    @debug padamts

    u0 = dict([k => zeros(FF, sz...) for k = all_field_names])
    # add_current_keys!(u0)
    u0 = groupkeys(u0)

    if npoles == 0
        u0[:Pn] = 0
        u0[:Jn] = 0
    else
        nE = length(Enames)
        u0[:Pn] = [[zeros(FF, sz...) for _ = 1:nE] for _ = 1:npoles]
        u0[:Jn] = [[zeros(FF, sz...) for _ = 1:nE] for _ = 1:npoles]
    end
    u0 = namedtuple(u0)

    centers = vmap(offsets) do v
        start = 1.5 + v |> FF
        stop = length.(rulers) - 0.5 + v |> FF
        n = round.(Int, stop - start) + 1
        getindexf.(rulers, range.(start, stop, n))
    end

    diffdeltas = kvmap(centers) do k, v
        k, map(v, rulers, offsets[k], 1:N) do c, r, o, i
            v = if o == 0
                b = boundaries[i, 1]
                if b == :periodic
                    x = r[end] - c[end]
                else
                    x = c[1] - r[1]
                end
                diff(vcat([r[1] - x], c))
            elseif o == -0.5
                b = boundaries[i, 2]
                if b == :periodic
                    x = c[1] - r[1]
                else
                    x = r[end] - c[end]
                end
                diff(vcat(c, [r[end] + x]))
            else
                error("Invalid offset")
            end
            s = (1:N) .== i
            reshape(v, (.!s + length(v) * s)...)
        end
    end
    # end

    sz = Tuple(sz)
    deltas = diff.(rulers)
    weights = map(prod, Base.product(deltas...))

    @debug offsets
    @debug padvals
    @debug padamts
    @debug boundvals

    geometry = pad_geometry(geometry, padvals, padamts) |> pairs |> OrderedDict
    centgeom = pad_geometry(centgeom, padvals, padamts) |> pairs |> OrderedDict
    centgeom = vmap(first, centgeom)

    # a = geometry[:σ][1]
    # fig = heatmap(a[:, :, size(a, 3)÷2])
    # CairoMakie.save("a.png", fig)
    # # GLMakie.volume(geometry[:σ][1],) |> display
    # error()

    for k = geometry_names
        @debug k extrema(geometry[k][1])
        if k == :invϵ
            @debug extrema(geometry[k][1]), extrema(geometry[k][2]), extrema(geometry[k][3])
        end
    end

    d = diffdeltas(r"H(.)") |> values
    b = boundvals(r"H(.)") |> values
    ∇ₕ = Del(
        map(1:N) do i
            getindex.(d, i) .|> array
        end,
        map(1:N) do i
            getindex.(b, i, 1)
        end,
        map(1:N) do i
            getindex.(b, i, 2)
        end,)

    d = diffdeltas(r"E(.)") |> values
    b = boundvals(r"E(.)") |> values
    ∇ₑ = Del(
        map(1:N) do i
            getindex.(d, i) .|> array
        end,
        map(1:N) do i
            getindex.(b, i, 1)
        end,
        map(1:N) do i
            getindex.(b, i, 2)
        end,)
    Ibbox1 = Ibbox = (:).(1 + padamts[:, 1], sz - padamts[:, 2])
    # @show Ibbox, padamts, sz
    # error("debug Ibbox")
    # Ibbox1 = (:).(max.(1, padamts[:, 1]), min.(sz, sz - padamts[:, 2] + 1))
    # Ibbox1 = (:).(first.(Ibbox) - nadj[:, 1], last.(Ibbox) + nadj[:, 2])

    grid = (; boundaries, sz, rulers, deltas, weights, bbox, L, offsets, searchers, boundvals, diffdeltas, padvals, ∇ₕ, ∇ₑ, padamts, Ibbox, Ibbox1, N, field_names, F, all_field_names, dt, bg, funnel) |> pairs |> OrderedDict
    writejsonnp(joinpath(path, "grid.json"), merge((; rulers, bbox, L) |> cpu |> FF, (; N, Ibbox, boundaries)))

    println("making sources...")
    mode_solutions = map(modes) do mode
        @unpack wavelengths, ports, voltage_line, current_loop, L = mode
        wavelengths, voltage_line, current_loop, L = FF.((wavelengths, voltage_line, current_loop, L))

        modes = mode(:modes)
        plane_size = size(modes[1][1](1))
        modes = map(modes) do v
                    map(v) do d
                        kvmap(d) do k, v
                            if all(==(v[1]), v)
                                v = v[1]
                            end
                            Symbol(k) => FF(v)
                        end
                    end
                end |> stack |> permutedims

        wavelengths /= λ
        L /= λ
        voltage_line /= λ
        current_loop /= λ

        frequencies = 1 ./ wavelengths
        dA = prod(L) / prod(plane_size)
        dr = L ./ plane_size

        # heatmap(real(modes[1].Ey)) |> display
        # error("stop")

        (; modes, wavelengths, frequencies, ports, voltage_line, current_loop, plane_size, dA, dr)
    end

    source_instances = SourceInstance.(sources, λ, (λs,), (grid,), (ϵ,); z, mode_solutions,)
    println("making monitors...")
    monitor_instances = MonitorInstance.(monitors, λ, (λs,), (grid,), (ϵ,); z, mode_solutions,)
    println("making designs...")
    if isempty(designs)
        design_instances = []
    else
        design_instances = DesignInstance.(designs, λ, (grid,), [centgeom];)
    end
    dimensions = last.(rulers) - first.(rulers)
    source_duration = maximum(duration.(source_instances))

    if isnothing(tmax)
        if !isnothing(path_length_multiple)
            tmax = path_length_multiple * sum(dimensions * nmax_dielectric) + source_duration
        else
            if isnothing(energy_decay_threshold)
                energy_decay_threshold = 0.03
            end
            tmax = 1000
        end
    end
    tmax += source_duration
    # @show tmax, dt
    tmax = dt * round(tmax / dt)

    nt = round(Int, FF(tmax) / dt)
    if !isempty(design_instances)
        gg = 1# max(1, round(Int, sqrt(nt / 100)))
        if isnothing(gradckptat)
            gradckptat = gg
        end
        nt = gradckptat * ceil(Int, nt / gradckptat)
        tmax = dt * nt
    end
    tmax = FF(tmax)

    nv = prod(sz)
    load = nt * nv

    λs /= λ
    λs = sort(λs)
    fs = (1 ./ λs)

    geometry = geometry |> 𝐟
    source_instances = fmap.(𝐟, source_instances)
    ∇ₕ = fmap(𝐟, ∇ₕ)
    ∇ₑ = fmap(𝐟, ∇ₑ)

    prob = (;
               grid, λ, λs, fs, relative_courant,
               source_instances, monitor_instances, design_instances,
               field_names, approx_2D_mode,
               tmax, energy_decay_threshold,
               geometry, centgeom, meshprops,
               u0, path, dt,
               saveat, checkat, gradckptat, verbosity, hasPEC,
               designs, targets, array, umax=OrderedDict(),
               name, time_extrapolation, timestamp,
           ) |> pairs |> OrderedDict
    if array == Array
        backend = :CPU
    else
        backend = :GPU
        for k in (:u0, :geometry, :source_instances, :monitor_instances, :fs)
            prob[k] = gpu(array, prob[k],)
        end
        for k = [:diffdeltas, :deltas, :weights]
            prob[:grid][k] = gpu(array, prob[:grid][k])
        end
    end
    println("""

    $BREAK
    simulation config

    backend: $backend
    float: $FF

    original size: $sz0
    padded size: $(sz)
    cell count: $(nv|>disp)

    relative Courant number: $relative_courant
    step size: $(1/dt|>round) steps/period
    source duration: $(source_duration|>disp) periods
    max time: $(tmax|>disp) periods    """)

    if isnothing(energy_decay_threshold)
        println("""
        time steps: $(nt|>disp)
        computation load: $(load|>disp) cell-steps
        """)
    else
        println("energy decay threshold: $(energy_decay_threshold)")
    end

    # if !pro && (nv > 1e6)
    #     println("Error: The free version of Lumi is limited to 1M cells for simulation and 10M cell-steps for inverse design. Consider upgrading to the pro version for larger simulations.\n")
    #     exit(2)
    # end
    prob
end
# update = update
# setup = setup

# dx = mean(cummin.(deltas, 4))
# v = min(0.8 * dx / dt, 0.1 / dt) |> FF
# lpml = -log(1.0f-2) / nmin / (2v) |> FF
# σpml = ϵmin * v
# mpml = μmin * v
# @debug (; dt, lpml, σpml, mpml, nmin)
# for i = 1:N
#     select = i .== 1:N
#     xyz = para = perp = [:x, :y, :z]
#     perp = [popat!(para, i)]
#     for j = 1:2
#         b = boundaries[i, j]
#         relpml = relative_pml_depths(i)(j)
#         if b == :PML
#             d = j == 1 ? deltas[i][1] : deltas[i][end]
#             npml = lpml / d
#             npml *= relpml
#             npml = max(ceil(Int, npml), 1)
#             start = !flush[i, j]
#             npml += flush[i, j]
#             v = d * (1:npml)
#             if j == 1
#                 # pushfirst!(deltas[i], fill(maxdeltas[i], npml)...)
#                 prepend!(rulers[i], rulers[i][1] - reverse(v))
#             else
#                 # push!(deltas[i], fill(maxdeltas[i], npml)...)
#                 append!(rulers[i], rulers[i][end] + v)
#             end

#             padamts(:default)[i, j] = npml
#             for k = geometry_names
#                 padamts[k][i, j] = npml
#             end

#             for k = (:ϵ, :μ, :invϵ, :invμ)
#                 padvals[k][i, j] = :replicate
#             end

#             padvals[:σ][i, j] = Ramp(σpml * relpml, start)
#             padvals[:m][i, j] = Ramp(mpml * relpml, start)

#             sz[i] += npml
#         end
#     end
# end
