using Luminescent, Random
Random.seed!(1)

push!(ARGS, lastrun("sparams"))
Luminescent.julia_main()
pop!(ARGS)
push!(ARGS, lastrun("inverse_design"))
Luminescent.julia_main()
