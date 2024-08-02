using PackageCompiler

create_app(".", "../LuminescentCompiled",
    executables=[
        "lumi" => "julia_main",
        #     "_1" => "julia_main1",
        #     "_2" => "julia_main2",
    ],
    precompile_execution_file="build/precompile_app.jl",
    include_transitive_dependencies=false,
    force=true
)