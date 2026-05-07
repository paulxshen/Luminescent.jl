

# function lumicu()::Cint
#     if !isempty(ARGS)
#         picrun(ARGS[1]; array=cu)
#     end
#     return 0
# end

function lastrun(name=nothing; study=nothing, wd=joinpath(pwd(), "studies"))
    if !isnothing(name)
        p = joinpath(wd, name)
        println("using simulation folder $p")
        return p
    end

    l = filter(isdir, readdir(wd, join=true))
    sort!(l, by=p -> Dates.unix2datetime(mtime(p)), rev=true)

    if !isnothing(study)
        for p = l
            try
                if open(joinpath(p, "solution.json")) do f
                    JSON.parse(f)["study"]
                end == study
                    println("using simulation folder $p")
                    return p
                end
            catch e
                println(e)
            end
        end
    end
    p = l[1]
    println("using simulation folder $p")
    return p
end
function parsekey(k)
    k = string(k)
    contains(k, ",") && return parsekey.(split(k, ","))
    if contains(k, "@")
        p, m = split(k, "@")
        return Symbol(p), parse(Int, m)
    end
    # parse(Int, k), 0
    Symbol("o$k"), 0
end
function getparam(sol, t, k, monitor_instances; f=nothing, λ=nothing)
    t = Symbol(t)
    if t == :S || t == :T
        s = sol.waves
    elseif t == :flux
        s = sol.profiles
    elseif t == :delay
        s = sol.delays
        k = @ignore_derivatives parse(Int, k)
        return s[k]
    else
        error("unknown parameter type $t")
    end

    (p2, m2), (p1, m1) = @ignore_derivatives parsekey(k)
    p1 = findfirst(m -> string(m.name) == string(p1), monitor_instances)
    p2 = findfirst(m -> string(m.name) == string(p2), monitor_instances)
    (isnothing(p1) || isnothing(p2)) && return

    if isnothing(f) && isnothing(λ)
        if t == :S
            return map(values(s[p2])) do v
                v[m2+1][1]
            end ./ map(values(s[p1])) do v
                v[m1+1][2]
            end
        elseif t == :flux
            return (getindex.(s[p2], 1)) ./ (getindex.(s[p1], 2))
            # return -s[p2] ./ s[p1]
        end
    end

    if isnothing(f)
        f = 1 / λ
    end
    i1 = @ignore_derivatives findfirst(isapprox(f), frequencies(monitor_instances[p1]))
    i2 = @ignore_derivatives findfirst(isapprox(f), frequencies(monitor_instances[p2]))
    S = s[p2][i2][m2+1][1] / s[p1][i1][m1+1][2]
    # S = -s[p2][i2][m2+1] / s[p1][i1][m1+1]
    t == :T && return abs2(S)
end
function make_geometry(masks, margins, lb, dl, geometry, designs, design_config, materials; ϵeff=nothing, perturb=nothing)
    isnothing(masks) && return geometry
    mat = design_config.fill.material
    namedtuple([k => begin
        a = geometry[k]
        T = typeof(a)
        if k == :ϵ
            f = @ignore_derivatives materials(mat)(:epsilon) |> FF
            v = @ignore_derivatives minimum(a)
            if perturb == k
                f *= convert.(FF, 1.001)
            end

            b = Zygote.Buffer(a)
            copyto!(b, a)

            for (mask, margin, design) in zip(masks, margins, designs)

                o = round.(Int, (design.bbox[1] - lb) / dl) + 1
                if ndims(b) == 3
                    o = [o..., 1 + round(Int, (zcore - zmin) / dl)]
                    mask = stack(fill(mask, round(Int, thickness / dl)))
                else
                    # f = ϵeff[2] |> FF
                end
                o -= margin
                mask = mask * (f - v) + v
                # b[range.(o, o .+ size(mask) .- 1)...] = mask
                b[[i:j for (i, j) = zip(o, o .+ size(mask) .- 1)]...] = mask
            end
            copy(b) |> constructor(T)
        else
            # println("no fill for $k")
            a
        end
    end for k = keys(geometry)])
end
# using GLMakie: volume
function vis(sol, prob, field=:Hz, path=".")
    @unpack u, p, _p = sol |> cpu
    prob = prob |> cpu
    @unpack monitor_instances, source_instances, λ, = prob
    @unpack deltas, spacings, bbox, dl = prob.grid
    u = u(field)
    N = ndims(u)
    # if N == 3
    #     volume(u) |> display
    #     heatmap(u[:, :, round(Int, size(u, 3) / 2)]) |> display
    # else
    #     heatmap(u) |> display
    # end
    # return
    g = _p(:ϵ)
    u = upsample(u, spacings)
    # g = imresize(g, size(u))
    bbox /= dl
    ratio = int(deltas[1] / dl)

    plt = quickie(u, g; dl, λ, monitor_instances, ratio, source_instances, bbox)
    display(plt)

    try
        CairoMakie.save(joinpath(path, "run.png"), plt,)
    catch e
        println("save plot failed")
        println(e)
    end
end

function query(sol, q, i=(:))
    q = string(q)
    t = @ignore_derivatives q[1]
    q = @ignore_derivatives q[2:end]
    q_parts = @ignore_derivatives split(q, ",")
    q2, q1 = q_parts[1], q_parts[2]

    if @ignore_derivatives q1[1] in "0123456789"
        q1 = "o$(q1)"
    end
    q1 = "$(q1)-"
    if occursin("@", q1)
        y1 = sol[:waves][q1][i]
    else
        y1 = sol[:flux][q1][i]
    end

    if @ignore_derivatives q2[1] in "0123456789"
        q2 = "o$(q2)"
    end
    q2 = "$(q2)+"
    if occursin("@", q2)
        y2 = sol[:waves][q2][i]
    else
        y2 = sol[:flux][q2][i]
    end

    if t == 'S'
        return y2 ./ y1
    elseif t == 'T'
        if occursin("@", q1)
            y1 = abs2.(y1)
        end
        if occursin("@", q2)
            y2 = abs2.(y2)
        end
        return y2 ./ y1
    end
end