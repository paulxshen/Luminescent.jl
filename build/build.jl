using PackageCompiler
# ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""

# create_app(".", "../LuminescentAI",
create_app("dummy", "../LuminescentAI",
    # executables=[
    #     "lumi" => "lumi",
    # ],
    # precompile_execution_file="build/precompile_app.jl",
    force=true,
    # incremental=true,
    # include_lazy_artifacts=true, 
    # include_transitive_dependencies=false
)

