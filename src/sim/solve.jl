function bell(t, dt, u=nothing)
    ignore_derivatives() do

        # fmap(u) do x
        #     isnan(x) && error("NaN detected in field")
        # end

        if floor(t + 0.001) > floor(t - dt + 0.001)
            ENV["autodiff"] == "0" && println("simulation period $t, took $(timepassed()) seconds")
        end
        if !haskey(ENV, "t0")
            ENV["t0"] = time()
        end
    end
end

function f1(((u,), p, (dt, field_diffdeltas, field_diffpadvals, source_instances,)), t)
    bell(t, dt, u)
    # @time u = update(u, p, t, dt, field_diffdeltas, field_diffpadvals, source_instances)
    u = update(u, p, t, dt, field_diffdeltas, field_diffpadvals, source_instances;)
    ((u,), p, (dt, field_diffdeltas, field_diffpadvals, source_instances,))
end

function f2(((u, um), p, (dt, field_diffdeltas, field_diffpadvals, source_instances,), (t0, T, monitor_instances)), t)
    bell(t, dt, u)
    # @time u = update(u, p, t, dt, field_diffdeltas, field_diffpadvals, source_instances;)
    u = update(u, p, t, dt, field_diffdeltas, field_diffpadvals, source_instances;)
    ks = @ignore_derivatives [keys(u.E)..., keys(u.H)...]
    um += [
        begin
            um = namedtuple(ks .=> field.((u,), ks, (m,)))
            [
                begin
                    c = dt / T * cispi(-2(t - t0) / λ)
                    @nograd c
                    c * um
                end for λ = wavelengths(m)
            ]
        end for m = monitor_instances]
    ((u, um), p, (dt, field_diffdeltas, field_diffpadvals, source_instances,), (t0, T, monitor_instances))
end

function solve(prob, ;
    save_memory=false, ulims=(-3, 3), framerate=nothing, showfield=:Hz, path="", #subpixel=true,
    kwargs...)
    @unpack approx_2D_mode, dt, u0, geometry, _geometry, source_instances, monitor_instances, Ttrans, Tss, ϵeff, array = prob
    @unpack plane_deltas, F, N, sz, deltas, field_diffdeltas, field_diffpadvals, field_lims, dl, spacings, geometry_padvals, geometry_padamts, _geometry_padamts = prob.grid

    @nograd plane_deltas, dt, u0, source_instances, monitor_instances, Ttrans, Tss, ϵeff, F, N, sz, deltas, field_diffdeltas, field_diffpadvals, field_lims, dl, spacings, geometry_padvals, geometry_padamts, _geometry_padamts

    p = geometry
    _p = _geometry

    TEMP = joinpath(path, "temp")
    FIELDS = joinpath(path, "temp", "fields")
    mkpath(FIELDS)
    mkpath(TEMP)

    if !isnothing(framerate)
        npzwrite(joinpath(TEMP, "g.npy"), prob._geometry.ϵ)
    end

    println("preprocessing geometry...")
    p = pad_geometry(p, geometry_padvals, geometry_padamts)
    p = apply_subpixel_averaging(p, field_lims)

    _p = pad_geometry(_p, geometry_padvals, _geometry_padamts)
    ks = filter(k -> k != (:ϵ), keys(_p))
    if !isempty(ks)
        p1 = namedtuple([k => downsamplefield(_p[k] |> cpu, field_lims, spacings) for k = ks])
        p = merge(p, p1)
    end

    mesheps = _p.ϵ |> cpu
    _pec = mesheps .>= (PECVAL - TOL)
    if !any(_pec)
        println("no PEC regions found in geometry")
        invϵ = tensorinv(mesheps, field_lims, spacings)
        @assert eltype(eltype(invϵ)) == F
    else
        # global mesheps = pecfy(mesheps, _pec, field_lims, spacings)
        ϵ = downsamplefield(mesheps, field_lims, spacings)
        invϵ = map(ϵ) do a
            1 ./ a
        end
    end
    invμ = 1
    p = merge(p, (; invϵ, invμ))
    p = fmap(array, p, AbstractArray{<:Number})

    # @ignore_derivatives @show typeof(invϵ)
    # return sum(invϵ) |> sum

    durations = [Ttrans, Tss]
    T = cumsum(durations)
    us0 = (u0,)
    init = (us0, p, (dt, field_diffdeltas, field_diffpadvals, source_instances,))

    ts = 0:dt:T[1]-F(0.001)
    @nograd ts

    @ignore_derivatives delete!(ENV, "t0")

    println("propagating transient fields...")
    if save_memory
        (u,), = adjoint_reduce(f1, ts, init, ulims)
    else
        (u,), = reduce(ts; init) do us, t

            ignore() do
                if !isnothing(framerate) && t > 0
                    if t % (1 / framerate) < dt
                        (u,), = us
                        a = u(showfield)
                        npzwrite(joinpath(FIELDS, "$t.npy"), a)

                        # CairoMakie.save(joinpath(FIELDS, "$t.png"), quickie(a, g; monitor_instances, source_instances, ulims),)
                        # quickie(a, g; monitor_instances, source_instances)
                    end
                end
            end
            f1(us, t;)
        end
    end

    if array == Array
        fmap(u) do x
            isnan(x) && error("NaN detected in field")
        end
    end

    ts = ts[end]+dt:dt:T[2]-F(0.001)
    init = ((u, 0), p, (dt, field_diffdeltas, field_diffpadvals, source_instances,), (T[2], durations[2], monitor_instances,))
    @nograd ts


    println("accumulating dft fields...")
    if save_memory
        (u, um), = adjoint_reduce(f2, ts, init, ulims)
    else
        (u, um), = reduce(f2, ts; init)
    end
    # return (u.E.Ex + u.H.Hz + u.E.Ey) .|> abs |> sum

    ignore_derivatives() do
        t0 = parse(Float64, ENV["t0"])
        println("simulation done in $(time() - t0) seconds (includes some JIT compilation time).")
    end

    # return sum(abs, um[1][1].Ex + um[1][1].Ey + um[1][1].Hz)
    ulims = 0
    # @assert all([all(!isnan, a) for a = u])

    # volume(cpu(prob.monitor_instances[2].(1)) + 0.001ϵ(1) |> cpu) |> display
    # extrema(cpu(prob.monitor_instances[2].λmodes(1)(1).Ey))
    # extrema(cpu(prob.monitor_instances[2].λmodes(1)(1).Hx))
    # volume(cpu(abs.(prob.source_instances[1].sigmodes(1)[2].Jy))) |> display
    # for i = 1:2
    #     for k = (:Ey, :Hx)
    #         prob.monitor_instances[i]._λmodes(1)[1](k) |> extrema |> println
    #     end
    # end
    # global a = um
    # conj(a[1][1].Ey) .* a[1][1].Hx |> sum |> println
    # conj(a[2][1].Ey) .* a[2][1].Hx |> sum |> println
    # error()

    @nograd monitor_instances
    v = map(um, monitor_instances) do um, m
        map(um, wavelengths(m)) do um, λ
            um = localframe(um, m)
            ap = am = nothing
            if !isnothing(m.λmodes)
                md = first.(plane_deltas)
                modes = m.λmodes[λ]
                _modes = m._λmodes[λ]
                @nograd modes, _modes, md
                ap = inner.(modes, (um,), (md,))
                am = inner.(_modes, (um,), (md,))
            end

            um, ap, am
        end
    end

    # extrema(prob.monitor_instances[1].λmodes(1)(1).Hx)
    um = [[v[1] for v = v] for v in v]
    ap = [[v[2] for v = v] for v in v]
    am = [[v[3] for v = v] for v in v]

    # nm = length(monitor_instances)
    # nλ = length(wavelengths(monitor_instances[1]))
    # um = [[v[j][i][1] for j = 1:nλ] for i = 1:nm]
    # ap = [[v[j][i][2] for j = 1:nλ] for i = 1:nm]

    # am = [[v[j][i][3] for j = 1:nλ] for i = 1:nm]

    return Solution(u, p, _p, ulims, um, ap, am)
end


struct Solution
    u
    p
    _p
    ulims
    um
    ap
    am
end
@functor Solution

function (s::Solution)(k, m, w=1, mn=0)
    @unpack u, ulims, um, ap, am = s
    if k == "a+"
        return s.ap[m][w][mn+1]
    elseif k == "a-"
        return s.am[m][w][mn+1]
    elseif k == "P_TE"
        return flux(um[m][w], :TE)
    elseif k == "P_TM"
        return flux(um[m][w], :TM)
    elseif k == "P"
        return flux(um[m][w])
    elseif k == "um"
        return um[m][w]
    end
end
# heatmap(___p.invϵ[1, 1])

function getsparams(s, iλ=1, imode=1)
    @unpack u, ulims, um, ap, am = s
    # nλ = length(ap[1])
    np = length(ap)
    # ifelse.(1:np .== 1, am[j], ap[j]) 
    map(1:np) do i
        ap[i][iλ][imode]
    end / am[1][iλ][imode]
end