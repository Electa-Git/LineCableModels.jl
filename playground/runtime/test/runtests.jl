using Test, TOML, UUIDs
using LineCableModelsRuntime
import DBInterface, SQLite, RequiredInterfaces
const RT = LineCableModelsRuntime
include("profiles.jl")
include("owned_commands.jl")
include("terminal_process.jl")
include("container_host.jl")
include("resource_journal.jl")
include("container_recovery.jl")
include("container_policy.jl")
include("container_image_layout.jl")
include("managed_agent.jl")
include("native_recovery.jl")
include("native_policy.jl")
include("container_scientific_driver.jl")
include("managed_terminal_driver.jl")
include("environment_fingerprint.jl")
include("catalogue.jl")
include("architecture.jl")
include("runtime_artifacts.jl")
include("sql_rows.jl")
include("worker_store.jl")
include("worker_inventory.jl")
include("broker_policy.jl")
include("tls_reader.jl")
include("assignments.jl")
include("agent_leases.jl")
include("terminal_resources.jl")
include("agent_terminals.jl")
include("lease_coordinator.jl")
include("terminal_coordinator.jl")
include("terminal_sockets.jl")
include("job_store.jl")
include("job_cancellation_store.jl")
include("scientific_coordinator.jl")
include("control_config.jl")
include("job_coordinator.jl")
include("agent_service.jl")
include("scientific_resources.jl")
include("startup_compilation.jl")
include("agent_science.jl")
include("job_cancellation.jl")
include("agent_jobs.jl")

struct MissingApplication <: AbstractApplication end
struct CountingApplication <: AbstractApplication
    count::Base.RefValue{Int}
end
RT.describe(::CountingApplication) =
    ApplicationDefinition("count", "Count", :workbench, "/workbenches/count")
RT.ui_command(app::CountingApplication, ::HostContext) = (app.count[] += 1; `false`)

@testset "identity and request boundary" begin
    key = repeat("secret", 8)
    proxy = ProxyIdentity("https://lcm.test", ["127.0.0.1"], key; administrators=["operator"])
    headers = ["X-LCM-Proxy-Key" => key, "X-LCM-Principal" => "alice"]
    @test authenticate(proxy, headers, "127.0.0.1").id == "alice"
    @test !authenticate(proxy, headers, "127.0.0.1").administrator
    @test authenticate(proxy, [headers[1], "X-LCM-Principal" => "operator"], "127.0.0.1").administrator
    @test_throws AccessDenied authenticate(proxy, headers, "192.0.2.1")
    @test_throws AccessDenied authenticate(proxy, ["X-LCM-Principal" => "alice"], "127.0.0.1")
    @test_throws AccessDenied authenticate(proxy, [headers; "x-lcm-principal" => "operator"], "127.0.0.1")
    @test_throws AccessDenied authenticate(proxy, [headers[1], "X-LCM-Principal" => "../alice"], "127.0.0.1")
    @test !occursin(key, repr(proxy))
    @test !occursin(string(proxy.key_digest), repr(proxy))
    @test_throws ArgumentError ProxyIdentity("http://lcm.test", ["127.0.0.1"], key)
    @test_throws ArgumentError ProxyIdentity("https://lcm.test/", ["127.0.0.1"], key)
    @test_throws ArgumentError ProxyIdentity("https://lcm.test", ["localhost"], key)
    @test_throws ArgumentError ProxyIdentity("https://lcm.test", ["127.0.0.1"], "short")
    @test_throws ArgumentError Principal("a\nb")
    valid = [headers; "Origin" => "https://lcm.test"; "X-LCM-Request" => "1"]
    @test authorize_request(proxy, valid, "127.0.0.1"; method="POST").id == "alice"
    @test authorize_request(proxy, valid, "127.0.0.1"; websocket=true).id == "alice"
    @test_throws AccessDenied authorize_request(proxy, headers, "127.0.0.1"; method="POST")
    @test_throws AccessDenied authorize_request(proxy, headers, "127.0.0.1"; websocket=true)
    @test_throws AccessDenied authorize_request(proxy, [headers; "Origin" => "https://evil.test"; "X-LCM-Request" => "1"], "127.0.0.1"; method="POST")
    @test_throws AccessDenied authorize_request(proxy, [headers; "Origin" => "https://lcm.test"], "127.0.0.1"; method="POST")
    @test_throws AccessDenied authorize_request(proxy, [valid; "origin" => "https://lcm.test"], "127.0.0.1"; websocket=true)
    local_policy = LocalIdentity("http://127.0.0.1:8080", Principal("developer"))
    @test authenticate(local_policy, Pair{String,String}[], "127.0.0.1").id == "developer"
    @test_throws AccessDenied authenticate(local_policy, headers, "127.0.0.1")
    @test_throws AccessDenied authenticate(local_policy, Pair{String,String}[], "192.0.2.1")
    @test_throws ArgumentError LocalIdentity("http://example.test", Principal("developer"))
end

function local_config(; overrides=Dict{String,Any}())
    merge(Dict{String,Any}(
        "schema_version" => 1, "enabled" => false,
        "gateway" => Dict("public_origin" => "http://127.0.0.1:8080"),
        "identity" => Dict("mode" => "local-development", "principal" => "alice"),
    ), overrides)
end

include("config_templates.jl")

@testset "strict configuration is inert" begin
    mktempdir() do dir
        path = joinpath(dir, "runtime.toml")
        function parseconfig(data)
            open(path, "w") do io
                TOML.print(io, data)
            end
            return read_config(path)
        end
        config = parseconfig(local_config())
        @test !config.enabled
        @test config.database == joinpath(dir, "state", "runtime.sqlite")
        @test config.scratch_root == joinpath(dir, "state", "runs")
        @test !isdir(joinpath(dir, "state"))
        @test config.limits.max_runs == 8
        @test_throws ArgumentError parseconfig(local_config(overrides=Dict("schema_version" => true)))
        @test_throws ArgumentError parseconfig(local_config(overrides=Dict("schema_version" => 2)))
        @test_throws ArgumentError parseconfig(local_config(overrides=Dict("enabled" => "yes")))
        @test_throws ArgumentError parseconfig(local_config(overrides=Dict("typo" => 1)))
        @test_throws ArgumentError parseconfig(local_config(overrides=Dict("gateway" =>
            Dict("public_origin"=>"http://127.0.0.1:8080", "listen_host"=>"0.0.0.0"))))
        @test_throws ArgumentError parseconfig(local_config(overrides=Dict("storage"=>Dict("scratch_root"=>"."))))
        @test_throws ArgumentError RunLimits(max_runs=0)
        @test_throws ArgumentError RunLimits(max_runs=true)
        @test_throws ArgumentError RunLimits(max_runs=1, max_runs_per_owner=2)
        @test_throws ArgumentError RunLimits(startup_seconds=Inf)
        @test_throws ArgumentError RunLimits(shutdown_seconds=0)
        keyfile = joinpath(dir, "proxy.key")
        write(keyfile, repeat("x", 40))
        chmod(keyfile, 0o600)
        proxied = Dict("schema_version"=>1, "gateway"=>Dict("public_origin"=>"https://lcm.test"),
            "identity"=>Dict("mode"=>"proxy", "proxy_key_file"=>"proxy.key"))
        @test parseconfig(proxied).identity isa ProxyIdentity
        chmod(keyfile, 0o644)
        @test_throws ArgumentError parseconfig(proxied)
    end
end

@testset "cheap, required application registration" begin
    role = RuntimeRequirement("parameters", ("line-parameters",))
    definition = ApplicationDefinition("study", "Cable study", :workbench,
        "/workbenches/study"; requirements=(role,))
    @test definition.requirements == (role,)
    @test_throws ArgumentError RuntimeRequirement("role", ())
    @test_throws ArgumentError RuntimeRequirement("role", ("duplicate", "duplicate"))
    @test_throws ArgumentError ApplicationDefinition("bad", "Bad", :other, "/path")
    @test_throws ArgumentError ApplicationDefinition("bad", "Bad", :workbench, "//evil.test")
    @test_throws ArgumentError ApplicationDefinition("bad", "Bad", :workbench, "/../secret")
    @test_throws ArgumentError ApplicationDefinition("bad", "Bad", :workbench, "/path"; requirements=(role, role))
    registry = ApplicationRegistry()
    count = Ref(0)
    @test register!(registry, CountingApplication(count)).id == "count"
    @test count[] == 0
    @test_throws ArgumentError register!(registry, CountingApplication(count))
    @test_throws ArgumentError register!(registry, MissingApplication())
    @test count[] == 0
end

@testset "durable ownership and transactional admission" begin
    mktempdir() do dir
        path = joinpath(dir, "runtime.sqlite")
        store = RuntimeStore(path)
        alice, bob = Principal("alice"), Principal("bob")
        application = ApplicationDefinition("study", "Study", :workbench, "/study")
        limits = RunLimits(max_runs=2, max_runs_per_owner=1)
        request_id = uuid4()
        try
            first = reserve_run!(store, alice, application; limits, request_id)
            @test first.owner == "alice"
            @test first.state == :reserved
            @test reserve_run!(store, alice, application; limits, request_id).id == first.id
            @test_throws CapacityUnavailable reserve_run!(store, alice, application; limits)
            second = reserve_run!(store, bob, application; limits)
            @test first.id != second.id
            @test_throws CapacityUnavailable reserve_run!(store, Principal("charlie"), application; limits)
            @test_throws AccessDenied get_run(store, bob, first.id)
            @test_throws AccessDenied transition_run!(store, bob, first.id, :starting)
            @test length(list_runs(store, alice)) == 1
            @test length(list_runs(store, Principal("operator"; administrator=true))) == 2
            @test_throws ArgumentError transition_run!(store, alice, first.id, :running)
            @test transition_run!(store, alice, first.id, :starting).state == :starting
            @test transition_run!(store, alice, first.id, :running).state == :running
            @test transition_run!(store, alice, first.id, :running).state == :running
            @test transition_run!(store, alice, first.id, :stopping).state == :stopping
            @test transition_run!(store, alice, first.id, :stopped).state == :stopped
            @test_throws ArgumentError transition_run!(store, alice, first.id, :running)
            @test reserve_run!(store, alice, application; limits).state == :reserved
            @test_throws ArgumentError transition_run!(store, bob, second.id, :failed; reason="private\ntrace")
        finally
            close(store)
        end
        reopened = RuntimeStore(path)
        try
            @test length(list_runs(reopened, alice)) == 2
            @test length(list_runs(reopened, bob)) == 1
        finally
            close(reopened)
        end
        unknown = joinpath(dir, "future.sqlite")
        db = SQLite.DB(unknown)
        RT.sql_rows(db, "PRAGMA user_version = 99")
        close(db)
        @test_throws ArgumentError RuntimeStore(unknown)
        other = joinpath(dir, "other.sqlite")
        db = SQLite.DB(other)
        RT.sql_rows(db, "CREATE TABLE unrelated (id INTEGER)")
        close(db)
        @test_throws ArgumentError RuntimeStore(other)
    end
end

@testset "admission across connections" begin
    mktempdir() do dir
        path = joinpath(dir, "runtime.sqlite")
        stores = (RuntimeStore(path), RuntimeStore(path))
        application = ApplicationDefinition("study", "Study", :workbench, "/study")
        limits = RunLimits(max_runs=1, max_runs_per_owner=1)
        results = Channel{Any}(2)
        try
            @sync for (index, store) in enumerate(stores)
                Threads.@spawn try
                    put!(results, reserve_run!(store, Principal("owner$index"), application; limits))
                catch error
                    put!(results, error)
                end
            end
            outcomes = [take!(results), take!(results)]
            @test count(x -> x isa RunRecord, outcomes) == 1
            @test count(x -> x isa CapacityUnavailable, outcomes) == 1
        finally
            foreach(close, stores)
        end
    end
end
