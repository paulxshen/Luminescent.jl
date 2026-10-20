tuts = ["photonics/bend", "photonics/ring", "photonics/TFLN_coupler", "photonics/GC",
    "RF/FSS", "RF/stub_filter", "inverse_design/GC_topopt", "inverse_design/GC_invd"]
for t = tuts
    x = basename(t)
    s = read("../moon/tutorials/$t/$x.ipynb", String)
    s = replace(s, r"/home/[^/\"]+/moon/tutorials/[^/\"]+/[^/\"]+/" => "", r"/home/[^/\"]+/" => "~/")
    write("$x.ipynb", s)
end
