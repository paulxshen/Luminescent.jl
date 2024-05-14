# a = Ex |> vec
# a = mode.Ex |> vec
# f = lines(abs.(a))

a = [[1 2]]
f(a) = sum(permutedims.(a, ((2, 1),))[1])
gradiente = gradient(f, a)