using TOML
port = parse(Int, ARGS[1])
directory = abspath(ARGS[2])
isdir(directory) || error("Expected owned test directory")
configuration = Dict(
    "schema_version"=>1, "enabled"=>true,
    "gateway"=>Dict("listen_host"=>"127.0.0.1", "port"=>port, "public_origin"=>"http://127.0.0.1:$port"),
    "identity"=>Dict("mode"=>"local-development", "principal"=>"browser-fixture"),
    "publisher"=>Dict("site_directory"=>normpath(joinpath(@__DIR__, "..", "..", "_site"))),
    "storage"=>Dict("database"=>joinpath(directory, "runtime.sqlite"), "scratch_root"=>joinpath(directory, "hosts")),
    "limits"=>Dict("max_runs"=>3, "max_runs_per_owner"=>3, "startup_seconds"=>120,
        "shutdown_seconds"=>1, "disconnect_grace_seconds"=>60),
)
open(joinpath(directory, "runtime.toml"), "w") do io
    TOML.print(io, configuration)
end
