tuts = ["photonic_circuits/bend", "photonic_circuits/SOI_coupler", "photonic_circuits/ring",
    "photonic_circuits/TFLN_coupler", "photonic_circuits/GC",
    "metasurfaces/FSS", "metasurfaces/metalens", "RF_circuits/stub_filter", "inverse_design/GC_invd"]
hub = "https://colab.luminescentai.com/hub/user-redirect/lab/tree/tutorials"
for t = tuts
    x = basename(t)
    s = read("../moon/tutorials/$t/$x.ipynb", String)
    s = replace(s, r"/home/[^/\"]+/moon/tutorials/[^/\"]+/[^/\"]+/" => "", r"/home/[^/\"]+/" => "~/")
    # Lumi Colab excludes inverse design
    if !startswith(t, "inverse_design")
        # after the title cell so MyST still takes its heading as the page title
        link = """{"cell_type": "markdown", "metadata": {}, "source": ["[Run in Lumi Colab]($hub/$t/$x.ipynb)"]},"""
        s = replace(s, r"\"cells\":\s*\[\n( *)\{.*?\n\1\},"s => SubstitutionString("\\0" * link); count=1)
    end
    write("$x.ipynb", s)
end

# py api from docstrings
api = """
import inspect, luminescent as lumi
groups = {
    "Workflow": ["mesh", "setup", "solve", "load", "movie", "view"],
    "Materials": ["material", "PEC_material", "placeholder_material"],
    "Sources & monitors": ["source", "monitor", "gf_source", "gf_monitor", "gaussian_mode"],
    "Frames": ["plane_frame", "points_frame", "box_frame", "cage_frame", "spherical_frame"],
    "Inverse design": ["design", "target", "optimizer", "apply_design", "design_from_gds"],
    "Utilities": ["straight", "remote_solve", "load_sparams", "version"],
}
# internal params hidden from signatures
hide = {"setup": ["entries", "keys", "verbose", "design_ports", "pixel_size", "info"],
        "view": ["material_color_intensity", "material_color_map", "materials"]}
s = "# Python API\\n"
for g, names in groups.items():
    s += f"## {g}\\n"
    for n in names:
        f = getattr(lumi, n)
        sig = inspect.signature(f)
        sig = sig.replace(parameters=[p for p in sig.parameters.values() if p.name not in hide.get(n, [])])
        s += f"### {n}\\n```python\\nlumi.{n}{sig}\\n```\\n{inspect.cleandoc(f.__doc__ or '')}\\n\\n"
open("api.md", "w").write(s)
"""
run(`python3 -c $api`)
