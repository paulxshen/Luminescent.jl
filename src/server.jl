using Genie, Genie.Renderer.Json, Genie.Requests
using HTTP, UnPack
include("main.jl")

function run_fdtd_server()
    Genie.Configuration.config!(
        server_port=8975,
        server_host="127.0.0.1",
    )

    route("/", method=GET) do
        "Hello, World!"
    end

    route("/local", method=POST) do
        @unpack action, path = jsonpayload()
        action = string(action)
        if action == "solve"
            picrun(path)
        elseif action == "query"
        end
    end

    up(async=false)
end

run_fdtd_server()