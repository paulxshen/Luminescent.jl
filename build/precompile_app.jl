using Luminescent, Random, JSON
Random.seed!(1)

for p = readdir("build/precompile_execution", join=true)
    # for p = ["tiny", "tiny3", "back"]
    #     p = joinpath("runs", p)
    if !contains(string(p), "CUDA")
        push!(ARGS, p)
        Luminescent.julia_main()
        pop!(ARGS)
    end
end
