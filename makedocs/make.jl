using Documenter
# include("../src/Luminescent.jl")
# using .Luminescent
# ex = (
#     "periodic_scattering",
#     "quarter_wavelength_antenna",
#     "inverse_design_waveguide_bend")
# for fn = ex
#     s = read("$(@__DIR__)/../examples/$fn/$fn.jl", String) |> strip
#     s = replace(s, "\r\n" => "\n")
#     if endswith(s, "=#")
#         s = s[1:end-2]
#     else
#         s = s * "\n```"
#     end
#     startswith(s, "#=") && (s = s[3:end])
#     s = "# " * (replace(fn, "_" => " ") |> titlecase) * "\nComplete file at [examples folder](https://github.com/paulxshen/Luminescent.jl/tree/master/examples)\n\n$s"
#     # replace!(s, "\"\"\"\n\n" => "\njulia```\n", "\n\n\"\"\"" => "\n```\n")
#     # replace!(s, "#=" => "julia```", "=#" => "```")
#     # s="\n hide\n"
#     # s = replace(s, r"\n[^\n]+# hide\n" => "\n")
#     # s="\n hide\n"
#     s = replace(s, r".+(hide).+\n" => "\n")
#     s = replace(s, "#=" => "```", "=#" => "```julia")
#     s = replace(s, r"```julia(\n)+" => "```julia\n", r"(\n)+```" => "\n```")

#     open("$(@__DIR__)/src/$fn.md", "w") do f
#         write(f, s,)
#     end
# end
makedocs(
    sitename="Luminescent.jl",
    format=Documenter.HTML(),
    # modules=[Luminescent],
    pages=[
        "index.md",
        # "guide.md",
        # "Tutorials" => ex .* [".md"],
        # "people.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#