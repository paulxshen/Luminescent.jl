@load "$(@__DIR__)/model" model0 model
p0 = make_geometry(model0, μ, σ, σm)
p = make_geometry(model, μ, σ, σm)
runsavesim(model0, 16, "pre_training")
runsavesim(model, 16, "post_training")

f = Figure()
a, = volume(f[1, 1], ones(2, 2, 2))
rotate_cam!(a.scene, 45, 5, 5)