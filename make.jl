tuts = ["photonics/bend", "photonics/ring", "photonics/TFLN_coupler", "photonics/GC",
    "RF/FSS", "RF/stub_filter", "inverse_design/GC_invd"]
hub = "https://hub.luminescentai.com/hub/user-redirect/lab/tree/tutorials"
for t = tuts
    x = basename(t)
    s = read("../moon/tutorials/$t/$x.ipynb", String)
    s = replace(s, r"/home/[^/\"]+/moon/tutorials/[^/\"]+/[^/\"]+/" => "", r"/home/[^/\"]+/" => "~/")
    # Lumi Hub excludes inverse design
    if !startswith(t, "inverse_design")
        link = """{"cell_type": "markdown", "metadata": {}, "source": ["[Run in Lumi Hub]($hub/$t/$x.ipynb)"]},"""
        s = replace(s, r"\"cells\":\s*\[" => SubstitutionString("\\0" * link); count=1)
    end
    write("$x.ipynb", s)
end
