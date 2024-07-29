println("starting optimization... first iter will be slow due to compilation.")
stop = false
img = nothing
for i = 1:maxiters
    global img = if virgin
        global virgin = false
        "before.png"
    elseif i == maxiters || stop
        "after.png"
    else
        img
    end
    @time global l, (dldm,) = withgradient(model) do m
        global sparams = write_sparams(m; img, autodiff=true)
        params = if target_type == "sparams"
            sparams
        else
            Porcupine.apply(abs2, sparams)
        end
        loss(params)
        # abs(sparams[1].mode_coeffs[2][1][1][1])

    end
    if stop
        break
    end
    if i == 1
        global best0 = best = l
        global sparams0 = deepcopy(sparams)
    end
    # l < 0 && break
    Flux.update!(opt_state, model, dldm)
    if l < best
        global best = l
    end
    if l < minloss
        # @info "Loss below threshold, stopping optimization."
        println("Loss below threshold, stopping optimization.")
        global stop = true
        global img = "after.png"
    end
    if i % 15 == 0
        if best - best0 > -0.01
            println("Loss stagnating, stopping optimization.")
            global stop = true
        else
            global best0 = best
        end
    end
    println("$i loss $l\n")
end
# @info "Done in $(time() - t0) ."
println("Done in $(time() - t0) .")
for (i, (m, d)) = enumerate(zip(model, designs))
    Images.save(joinpath(path, "design$i.png"), Gray.(m() .< 0.5))
end
sol = (;
    before=sparam_family(sparams0),
    after=sparam_family(sparams),
    designs, path, dx,
    designs=[m() .> 0.5 for m in model],
)
# @save "$path/sol.json" sol
open("$(path)/sol.json", "w") do f
    write(f, json(sol))
end