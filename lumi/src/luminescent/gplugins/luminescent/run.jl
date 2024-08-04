
using Luminescent, Random

Random.seed!(1)
if isempty(ARGS)
    path = lastrun()
    println("path: ", path)
else
    path = ARGS[1]
end
# gfrun(path; Courant=0.5)
gfrun(path)

# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl"
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl"
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl"