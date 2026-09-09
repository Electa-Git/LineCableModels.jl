using Test, UUIDs
using LineCableModelsRuntime
const RecoveryRT = LineCableModelsRuntime

function recovery_fence()
    RecoveryRT.Protocol.AssignmentFence(string(uuid4()), string(uuid4()), "alice", "main",
        "worker-a", string(uuid4()), string(uuid4()), "line-parameters", "1.0.0", repeat("a",64), 1)
end

function with_recovery_fixture(action::Function, kind::Symbol)
    mktempdir() do directory
        runner = CommandRunner()
        journal = ResourceJournal(joinpath(directory, "journal"), "worker-a")
        engine = ContainerEngine(kind, string(kind), true)
        info = kind == :podman ? Dict{String,Any}(
            "host"=>Dict("os"=>"linux", "serviceIsRemote"=>false, "cgroupVersion"=>"v2",
                "cgroupControllers"=>["memory", "pids"],
                "security"=>Dict("seccompEnabled"=>true, "rootless"=>true)),
            "store"=>Dict("graphRoot"=>"/owned/storage")) : Dict{String,Any}(
            "ID"=>"daemon-one", "OSType"=>"linux", "CgroupVersion"=>"2",
            "MemoryLimit"=>true, "SwapLimit"=>true, "CpuCfsPeriod"=>true,
            "CpuCfsQuota"=>true, "PidsLimit"=>true, "SecurityOptions"=>["name=seccomp,profile=builtin"])
        object = RecoveryRT.JSON3.read(RecoveryRT.JSON3.write(info))
        host = kind == :podman ? RecoveryRT.podman_host_check(engine, object) :
            RecoveryRT.docker_host_check(engine, object, "unix:///owned/docker.sock")
        state = Dict{Symbol,Any}(:inspect=>true, :list=>true, :remove=>true, :malformed_list=>false,
            :change_after_stop=>false, :storage_inode=>"42", :scope_reads=>0)
        probe = args -> begin
            state[:scope_reads] += 1
            (true, RecoveryRT.JSON3.write(info))
        end
        machine_id = repeat("1",32)
        storage_identity = path -> (path, "1", state[:storage_inode])
        options = (; probe, machine_id, storage_identity)
        scope = container_scope(runner, host; options...)
        containers = Dict{String,Any}()
        history = Vector{String}[]
        change_scope!() = kind == :podman ? (state[:storage_inode] = "43") : (info["ID"] = "daemon-two")
        restore_scope!() = kind == :podman ? (state[:storage_inode] = "42") : (info["ID"] = "daemon-one")
        function invoke(args)
            push!(history, copy(args))
            offset = findfirst(==("container"), args)
            offset === nothing && error("unexpected fixture command")
            action = args[offset+1]
            id = last(args)
            if action == "inspect"
                state[:inspect] || return CommandResult(1, "", "private-inspect-failure")
                found = findfirst(pair -> first(pair) == id || last(pair)["Name"] == id, collect(containers))
                found === nothing && return CommandResult(1, "", "private-not-found")
                value = collect(containers)[found].second
                return CommandResult(0, RecoveryRT.JSON3.write(value), "")
            elseif action == "ls"
                state[:list] || return CommandResult(1, "", "private-list-failure")
                state[:malformed_list] && return CommandResult(0, "invalid private inventory", "")
                filter = args[findfirst(==("--filter"), args)+1]
                lines = String[]
                for (key, value) in containers
                    if startswith(filter, "id=") ? key == filter[4:end] : value["Name"] == filter[6:end]
                        push!(lines, key * " " * value["Name"])
                    end
                end
                return CommandResult(0, join(lines, "\n"), "")
            elseif action == "stop"
                state[:change_after_stop] && change_scope!()
                return CommandResult(0, "", "")
            elseif action == "rm"
                state[:remove] || return CommandResult(1, "", "private-remove-failure")
                pop!(containers, id, nothing)
                return CommandResult(0, "", "")
            end
            error("unexpected fixture command")
        end
        function add!(; bound=false)
            receipt = reserve_resource!(journal, recovery_fence(), kind, scope)
            id = bytes2hex(RecoveryRT.sha256(string(uuid4())))
            containers[id] = Dict("Id"=>id, "Name"=>resource_name(receipt),
                "Config"=>Dict("Labels"=>resource_labels(receipt)))
            return bound ? bind_resource!(journal, receipt, id) : receipt, id
        end
        try
            action((; runner, journal, host, info, scope, options=(;options..., invoke),
                containers, history, state, add!, change_scope!, restore_scope!))
        finally
            close(runner)
            close(journal)
        end
    end
end

@testset "both engine adapters retire exact owned IDs and verify absence" begin
    for kind in (:podman, :docker)
        with_recovery_fixture(kind) do f
            (; receipt, id) = let (receipt, id) = f.add!(); (; receipt, id) end
            other, unrelated = f.add!(; bound=true)
            # The first resource is a crash-gap intent; cleanup must bind its full ID.
            @test remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test !haskey(f.containers, id)
            @test haskey(f.containers, unrelated)
            @test only(resource_receipts(f.journal)).id == other.id
            mutations = filter(args -> "stop" in args || "rm" in args, f.history)
            @test length(mutations) == 2
            @test all(args -> last(args) == id, mutations)
            @test !any(args -> any(in(("prune", "--all", "--volumes")), args), mutations)
            @test remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test length(filter(args -> "rm" in args, f.history)) == 1
            @test recover_containers!(f.journal, f.runner, f.host; f.options...) === nothing
            @test isempty(f.containers)
            @test isempty(resource_receipts(f.journal))
            kind == :podman && @test :cpu_controller_missing in f.host.failures
        end
    end
end

@testset "failed inspection or removal never masquerades as cleanup" begin
    for kind in (:podman, :docker)
        with_recovery_fixture(kind) do f
            receipt, id = f.add!(; bound=true)
            f.state[:inspect] = false
            @test !remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test length(resource_receipts(f.journal)) == 1
            @test !any(args -> "rm" in args, f.history)
            f.state[:list] = false
            @test !remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            f.state[:list] = true
            f.state[:malformed_list] = true
            @test !remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            f.state[:malformed_list] = false
            f.state[:inspect] = true
            f.state[:remove] = false
            @test !remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test only(resource_receipts(f.journal)).physical_id == id
            f.state[:remove] = true
            @test remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test isempty(resource_receipts(f.journal))
        end
    end
end

@testset "engine replacement and foreign labels retain cleanup ownership" begin
    for kind in (:podman, :docker)
        with_recovery_fixture(kind) do f
            receipt, id = f.add!()
            f.change_scope!()
            @test_throws ArgumentError remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test isempty(f.history)
            @test length(resource_receipts(f.journal)) == 1
            f.restore_scope!()
            f.containers[id]["Config"]["Labels"]["org.linecablemodels.journal"] = string(uuid4())
            @test_throws ArgumentError remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test !any(args -> "stop" in args || "rm" in args, f.history)
            f.containers[id]["Config"]["Labels"] = resource_labels(receipt)
            f.state[:change_after_stop] = true
            @test_throws ArgumentError remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
            @test !any(args -> "rm" in args, f.history)
            @test only(resource_receipts(f.journal)).physical_id == id
            f.restore_scope!()
            f.state[:change_after_stop] = false
            @test remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
        end
    end
end

@testset "one unresolved receipt does not prevent independent recovery attempts" begin
    with_recovery_fixture(:podman) do f
        invalid, bad_id = f.add!()
        valid, good_id = f.add!()
        f.containers[bad_id]["Config"]["Labels"] = Dict()
        @test_throws ArgumentError recover_containers!(f.journal, f.runner, f.host; f.options...)
        @test haskey(f.containers, bad_id)
        @test !haskey(f.containers, good_id)
        @test only(resource_receipts(f.journal)).id == invalid.id
        @test !any(args -> ("stop" in args || "rm" in args) && last(args) == bad_id, f.history)
    end
end

@testset "recovery cannot delete an unrecorded resource or revive a lease" begin
    with_recovery_fixture(:docker) do f
        receipt, id = f.add!(; bound=true)
        object = deepcopy(f.containers[id])
        @test remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
        f.containers[id] = object # Simulate external restoration after its receipt was removed.
        before = length(f.history)
        @test !remove_owned_container!(f.journal, f.runner, f.host, receipt; f.options...)
        @test haskey(f.containers, id)
        @test !any(args -> "rm" in args || "stop" in args, f.history[before+1:end])
        @test isempty(resource_receipts(f.journal))
    end
end
