using Luminescent, Random
Random.seed!(1)

for p = filter(isdir, readdir("build/precompile_execution", join=true))
    push!(ARGS, p)
    Luminescent.julia_main()
    pop!(ARGS)
end
