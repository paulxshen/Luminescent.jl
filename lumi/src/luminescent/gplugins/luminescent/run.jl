using Luminescent, Random
Random.seed!(1)
if isempty(ARGS)
    path = lastrun()
    println("path: ", path)
else
    path = ARGS[1]
end
sol = gfrun(path)
# @show sol.tparams

# using Pkg
# pkg"add Porcupine,Jello,ArrayPadding;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl;dev C:\Users\pxshe\OneDrive\Desktop\Luminescent.jl;up"
