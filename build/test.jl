# using Luminescent
include("../src/main.jl")
using Random
Random.seed!(1)

l = readdir("build/precompile_execution", join=true)
sort!(l, by=p -> Dates.unix2datetime(mtime(p)))

for path = filter(isdir, l)
    sol = gfrun(path)
end
