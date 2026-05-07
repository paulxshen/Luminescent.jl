include("C:\\Users\\default.LAPTOP-HMRU58MH\\Desktop\\lumi\\src\\main.jl")

studies = "C:\\Users\\default.LAPTOP-HMRU58MH\\Desktop\\lumi\\luminescent\\studies"
ENV["JULIA_DEBUG"] = "Main"


# for x = ["photonic_waveguide_bend", "microstrip_stub_filter", "microwave_frequency_selective_surface", "microstrip_patch_antenna"]
#     x = joinpath(studies, x)
#     picrun(x; gpu_backend="CUDA")
#     makemovie(x)
# end

picrun(joinpath(studies, "mux"))
