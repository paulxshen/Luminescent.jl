function print(x...)
    Base.print(x...)
    (haskey(ENV, "server") && ENV["server"] != "local") && @ignore_derivatives log_write(join(x, ""))
end
println(x...) = print(x..., "\n")
function error(x...)
    println(x...)
    Base.error(x...)
end

function _versioninfo()
    println("""
    CPU:
    \tthreads: $(Threads.nthreads())
    """)
    # versioninfo()
    # println()
    # Check if CUDA is available
    # if CUDA.functional()
    #     # Print general CUDA version information
    #     # CUDA.versioninfo()

    #     # Iterate through available CUDA devices and print their names
    #     for i in 0:CUDA.ndevices()-1
    #         device = CuDevice(i)
    #         println("""
    #         GPU $i: $(CUDA.name(device))
    #         \tCUDA Compute Capability: $(CUDA.capability(device))
    #         \tVRAM: $(CUDA.totalmem(device) / (1024^3)) GB
    #         """)
    #     end
    # else
    #     println("CUDA is not functional. Please ensure drivers and toolkit are installed correctly.")
    # end
end



# Function: log_write - add texts as rows in texts table
function log_write(text::String)
    job = ENV["job"]
    cloud = ENV["cloud"]
    # db = DBInterface.connect(MySQL.Connection, "136.116.160.55", "root", "balmy", db=ENV["cloud"],)
    # DBInterface.execute(db, "INSERT INTO log$job (val) VALUES ('$text')")
    run(`python3 $(ENV["pysql"]) log_write $cloud $job $text`)
end

function check_ping()
    if haskey(ENV, "server") && ENV["server"] != "local"
        ping = time()
        job = ENV["job"]
        clstr = ENV["cloud"]

        # db = DBInterface.connect(MySQL.Connection, "136.116.160.55", "root", "balmy", db=clstr)

        # r = DBInterface.execute(db, "SELECT val FROM state$job WHERE var = 'client_ping'") |> DataFrame
        # client_ping = parse(Float64, r[1, :val])

        # r = DBInterface.execute(db, "SELECT val FROM state$job WHERE var = 'status'") |> DataFrame
        # status = r[1, :val]

        # r = DBInterface.execute(db, "SELECT val FROM state$job WHERE var = 'timestamp'") |> DataFrame
        # timestamp = r[1, :val]
        # close(db)
        client_ping, status, timestamp = read(`python3 $(ENV["pysql"]) get_job_state $clstr $job`, String) |> x -> split(x, ",") |> x -> (parse(Float64, x[1]), x[2], x[3])

        (ping - client_ping > 30) && error("client disconnected")
        ping - parse(Float64, ENV["start"]) > 900 && error("exceeded max run time of 15 minutes")
        # status != "running" && error("job status $status. job aborted")
        ENV["timestamp"] != string(timestamp) && error("timestamp mismatch. newer job of same study detected. this job aborted")

        @debug "ping check passed"
    end
end
# check_ping() = retry(_check_ping; delays)()

function set_job_status(status)
    job = ENV["job"]
    cloud = ENV["cloud"]
    # db = DBInterface.connect(MySQL.Connection, "136.116.160.55", "root", "balmy", db=clstr,)
    # DBInterface.execute(db, "UPDATE state$job SET val = '$s' WHERE var = 'status'")
    run(`python3 $(ENV["pysql"]) set_job_status $cloud $job $status`)
end
function set_server_status(status)
    server = ENV["server"]
    # db = DBInterface.connect(MySQL.Connection, "136.116.160.55", "root", "balmy", db="servers",)
    # DBInterface.execute(db, "UPDATE $server SET val = '$s' WHERE var = 'status'")
    run(`python3 $(ENV["pysql"]) set_server_status servers $server $status`)
end