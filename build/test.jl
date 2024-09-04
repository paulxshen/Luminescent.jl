using Luminescent, Random
Random.seed!(1)

for path = filter(isdir, readdir("build/precompile_execution", join=true))
    sol = gfrun(path)
end
