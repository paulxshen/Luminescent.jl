using Optimisers
using SparseArrays
using Flux: Adam
for T in (:Float32, :Float16, :Float64)
    @eval Optimisers.maywrite(::CUDA.CUSPARSE.CuSparseMatrixCSC{$T,Int32}) = true
    @eval Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{$T,Int32}) = true
end

function julia_main()::Cint
    if !isempty(ARGS)
        gfrun(ARGS[1])
    end
    return 0
end

function lastrun(; name=nothing, study=nothing, wd="runs")
    path = joinpath(pwd(), wd)
    !isnothing(name) && return joinpath(path, name)

    l = filter(isdir, readdir(path, join=true))
    sort!(l, by=p -> Dates.unix2datetime(mtime(p)), rev=true)

    if !isnothing(study)
        for p = l
            try
                open(joinpath(p, "solution.json")) do f
                    JSON.parse(f)["study"]
                end == study && return p
            catch e
                println(e)
            end
        end
    end
    return l[1]
end

function write_sparams(runs, run_probs, origin, deltas,
    designs=nothing, design_config=nothing, models=nothing;
    alg=nothing, save_memory=false, verbose=false, perturb=nothing, framerate=0, path="", kw...)
    F = run_probs[1].F
    dx = deltas[1]
    sols = [
        begin
            prob[:_geometry] = make_geometry(models, origin, dx, prob._geometry, designs, design_config; F, perturb)
            #@debug typeof(prob.u0.E.Ex), typeof(prob.geometry.ϵ)
            sol = solve(prob; alg, save_memory, verbose, framerate, path)
        end for (i, prob) in enumerate(run_probs)
        # end for (i, prob) in enumerate(run_probs)
    ]
    # S = sols[1]("a+", 1) |> abs2
    # return (; S, sols)

    ulims = sols[1].ulims
    # return sol
    coeffs = OrderedDict()
    for (sol, run) in zip(sols, runs)
        sources = values(run.sources)
        monitors = values(run.monitors)
        source_port = first(sources).port
        # source_mn = first(sources).wavelength_mode_numbers(1)[1]
        source_mn = first(sources).wavelength_mode_numbers |> Porcupine.first |> Porcupine.first
        for (m, monitor) = enumerate(monitors)
            for (w, λ) = enumerate(keys(monitor.wavelength_mode_numbers))
                for mn = monitor.wavelength_mode_numbers[λ]
                    monitor_port = monitor.port
                    λ = Symbol(λ)
                    if !haskey(coeffs, λ)
                        coeffs[λ] = OrderedDict()
                    end
                    s = "$monitor_port@$mn," * "$source_port@$source_mn"
                    s = Symbol(s)
                    coeffs[λ][s] = (sol("a+", m, w, mn), sol("a-", m, w, mn))
                end
            end
        end
    end
    # return coeffs(1)(1)[1] |> abs2

    S = OrderedDict([λ => OrderedDict([k => begin
        s = ignore() do
            split(string(k), ",")[2]
        end
        # Symbol(
        coeffs[λ][k][1] / coeffs[λ][Symbol("$s,$s")][2]
    end for (k) = keys(coeffs[λ])]) for (λ) = keys(coeffs)])
    # if source_mn == mn == 0
    #     coeffs[λ]["$monitor_port,$source_port")] = v
    # end
    return (; S, sols)
end

function make_geometry(models, origin, dx, geometry, designs, design_config; F=Float32, perturb=nothing)
    isnothing(models) && return geometry
    ratio = models[1].ratio
    dx /= ratio
    namedtuple([k => begin
        if k in keys(design_config.fill)
            f = design_config.fill[k] |> F
            v = design_config.void[k] |> F
            if perturb == k
                f *= convert.(F, 1.001)
            end

            b = Zygote.Buffer(geometry[k])
            copyto!(b, geometry[k])

            for (m, design) in zip(models, designs)
                mask = m((x, r) -> x, v, f)
                # mask = m() * (f - v) + v

                o = round.(Int, (design.bbox[1] - origin) / dx) + 1
                if ndims(b) == 3
                    o = [o..., 1 + round(Int, (zcore - zmin) / dx)]
                    mask = stack(fill(mask, round(Int, thickness / dx)))
                end
                # b[range.(o, o .+ size(mask) .- 1)...] = mask
                b[[i:j for (i, j) = zip(o, o .+ size(mask) .- 1)]...] = mask
            end
            copy(b)
        else
            geometry[k]
        end
    end for k = keys(geometry)])
end
using GLMakie: volume
function plotsols(sols, probs, path, origin)
    for (i, (prob, sol)) in enumerate(zip(probs, sols))
        # try
        @unpack u, p, _p = sol |> cpu
        @unpack monitor_instances, source_instances, deltas, dl, λ, bbox = prob |> cpu
        u = u.Hz
        N = ndims(u)
        if N == 3
            volume(u) |> display
            heatmap(u[:, :, round(Int, size(u, 3) / 2)]) |> display
        else
            heatmap(u) |> display
        end
        return
        g = _p.ϵ
        u = imresize(u, Tuple(round.(Int, size(u) .* ratio)))
        g = imresize(g, size(u))
        bbox /= dl

        plt = quickie(u, g; dl, λ, monitor_instances, ratio, source_instances, origin, bbox)
        display(plt)

        if !isa(path, Base.AbstractVecOrTuple)
            path = (path,)
        end
        for path = path
            try
                CairoMakie.save(joinpath(path, "run_$i.png"), plt,)
            catch e
                println("save plot failed")
                println(e)
            end
        end
    end
end
