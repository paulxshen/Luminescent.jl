@show Threads.nthreads()
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
ENV["JULIA_DEBUG"] = "Main"

d = "~/moon/" |> expanduser
include("$d/src/main.jl")

name = "photonic_waveguide_bend"
name = "photonic_ring_resonator"
name = "microwave_FSS"
d = "~/moon/luminescent" |> expanduser
p = "$d/studies/$name"
picrun(p)
makemovie(p)
