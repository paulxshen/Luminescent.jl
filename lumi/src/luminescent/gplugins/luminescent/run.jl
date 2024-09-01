
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
# pkg"dev https://github.com/paulxshen/Luminescent.jl"
# pkg"add Porcupine,Jello,ArrayPadding;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl;dev C:\Users\pxshe\OneDrive\Desktop\Luminescent.jl;up"