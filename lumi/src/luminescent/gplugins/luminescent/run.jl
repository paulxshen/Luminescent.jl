
using Luminescent, Random
names(Luminescent)
Random.seed!(1)
if isempty(ARGS)
    path = lastrun()
    println("path: ", path)
else
    path = ARGS[1]
end
# gfrun(path; Courant=0.5)
sol = gfrun(path)
# @show sol.tparams
# # Pkg.resolve()

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl"
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl"
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl"
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Luminescent.jl"