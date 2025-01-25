using PackageCompiler
# ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""

# create_app(".", "../Luminescent",
#     # create_app("dummy", "../Luminescent",
#     # create_app("dummy", "../Luminescent550",
#     # executables=[
#     #     "lumi" => "lumi",
#     # ],
#     precompile_execution_file="build/precompile_app.jl",
#     force=true,
#     incremental=true,
#     # include_lazy_artifacts=true,
#     # include_transitive_dependencies=false
# )

create_sysimage(["Luminescent"]; precompile_execution_file="build/precompile_app.jl")