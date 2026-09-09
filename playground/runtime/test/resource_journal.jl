using Test, UUIDs
using LineCableModelsRuntime
const JournalRT = LineCableModelsRuntime
journal_fence(; worker="worker-a", lease=string(uuid4()), generation=1) =
    JournalRT.Protocol.AssignmentFence(lease, string(uuid4()), "alice", "main", worker,
        string(uuid4()), string(uuid4()), "line-parameters", "1.0.0", repeat("a",64), generation)
journal_changed(value::T; changes...) where T =
    T((get(changes, field, getfield(value, field)) for field in fieldnames(T))...)

@testset "resource receipts bind exact ownership without creating authority" begin
    mktempdir() do directory
        root = joinpath(directory, "journal")
        journal = ResourceJournal(root, "worker-a"; capacity=2)
        fence = journal_fence()
        try
            @test isempty(resource_receipts(journal))
            @test stat(root).mode & 0o777 == 0o700
            @test_throws ArgumentError ResourceJournal(root, "worker-a")
            @test_throws ArgumentError reserve_resource!(journal, journal_fence(worker="other"), :podman, repeat("e",64))
            @test_throws ArgumentError reserve_resource!(journal, fence, :unknown, repeat("e",64))
            intent = reserve_resource!(journal, fence, :podman, repeat("e",64))
            @test intent.physical_id === nothing
            @test intent.fence == fence
            @test only(resource_receipts(journal)).id == intent.id
            @test reserve_resource!(journal, fence, :podman, repeat("e",64)).id == intent.id
            @test_throws ArgumentError reserve_resource!(journal, fence, :docker, repeat("e",64))
            @test_throws ArgumentError reserve_resource!(journal, journal_changed(fence; generation=2), :podman, repeat("e",64))
            @test_throws ArgumentError reserve_resource!(journal, fence, :podman, repeat("f",64))
            @test_throws ArgumentError reserve_resource!(journal, fence, :podman, "unverified-scope")
            @test resource_name(intent) == "lcm-exec-" * string(intent.id)
            @test !occursin(root, repr(MIME"text/plain"(), journal))
            @test !occursin("alice", repr(MIME"text/plain"(), intent))
            @test_throws ArgumentError bind_resource!(journal, intent, "short-id")
            @test_throws ArgumentError bind_resource!(journal, intent, "../bad")
            bound = bind_resource!(journal, intent, repeat("b",64))
            @test bound.physical_id == repeat("b",64)
            @test bind_resource!(journal, intent, repeat("b",64)).physical_id == bound.physical_id
            @test_throws ArgumentError bind_resource!(journal, bound, repeat("c",64))
            @test_throws ArgumentError forget_resource!(journal, intent)
            labels = resource_labels(bound)
            @test matches_resource(bound, repeat("e",64), repeat("b",64), resource_name(bound), labels)
            @test !matches_resource(bound, repeat("f",64), repeat("b",64), resource_name(bound), labels)
            @test matches_resource(intent, repeat("e",64), repeat("b",64), resource_name(intent), labels)
            @test !matches_resource(bound, repeat("e",64), repeat("c",64), resource_name(bound), labels)
            @test !matches_resource(bound, repeat("e",64), repeat("b",64), "unrelated-container", labels)
            for key in keys(labels)
                changed = copy(labels)
                changed[key] = "unrelated"
                @test !matches_resource(bound, repeat("e",64), repeat("b",64), resource_name(bound), changed)
            end
            @test !matches_resource(bound, repeat("e",64), repeat("b",64), resource_name(bound), Dict())
            native = reserve_resource!(journal, journal_fence(), :native, repeat("e",64))
            @test endswith(resource_name(native), ".service")
            @test_throws ArgumentError bind_resource!(journal, native, repeat("b",64))
            native = bind_resource!(journal, native, repeat("c",32))
            @test native.physical_id == repeat("c",32)
            @test_throws ArgumentError reserve_resource!(journal, journal_fence(), :docker, repeat("e",64))
            @test forget_resource!(journal, bound) === nothing
            @test forget_resource!(journal, bound) === nothing
            @test length(resource_receipts(journal)) == 1
            @test forget_resource!(journal, native) === nothing
            @test isempty(resource_receipts(journal))
        finally
            close(journal)
        end
        @test close(journal) === nothing
        @test_throws ArgumentError resource_receipts(journal)
        @test_throws ArgumentError reserve_resource!(journal, fence, :podman, repeat("e",64))
        @test_throws ArgumentError ResourceJournal(root, "other-worker")
        restored = ResourceJournal(root, "worker-a")
        try
            @test restored.id == journal.id
            @test isempty(resource_receipts(restored))
        finally
            close(restored)
        end
    end
end

@testset "committed intent survives interrupted binding and stale writes are bounded" begin
    mktempdir() do directory
        root = joinpath(directory, "journal")
        journal = ResourceJournal(root, "worker-a")
        receipt = reserve_resource!(journal, journal_fence(), :docker, repeat("e",64))
        close(journal)
        pending = joinpath(root, ".pending-" * string(uuid4()))
        write(pending, "{incomplete atomic binding")
        chmod(pending, 0o600)
        restored = ResourceJournal(root, "worker-a")
        try
            @test !ispath(pending)
            recovered = only(resource_receipts(restored))
            @test recovered.id == receipt.id
            @test recovered.physical_id === nothing
            @test matches_resource(recovered, repeat("e",64), repeat("d",64), resource_name(recovered), resource_labels(recovered))
        finally
            close(restored)
        end
        @test isfile(joinpath(root, string(receipt.id) * ".json"))
    end
end

@testset "foreign, linked and malformed journal contents are never removed" begin
    mktempdir() do directory
        foreign = joinpath(directory, "foreign")
        mkdir(foreign; mode=0o700)
        write(joinpath(foreign, "user-file"), "untouched")
        @test_throws ArgumentError ResourceJournal(foreign, "worker-a")
        @test readdir(foreign) == ["user-file"]
        @test read(joinpath(foreign, "user-file"), String) == "untouched"
        broad = joinpath(directory, "broad")
        mkdir(broad; mode=0o755)
        # A caller's restrictive umask must not repair this rejection fixture.
        chmod(broad, 0o755)
        @test stat(broad).mode & 0o777 == 0o755
        @test_throws ArgumentError ResourceJournal(broad, "worker-a")
        link = joinpath(directory, "linked")
        symlink(foreign, link)
        @test_throws ArgumentError ResourceJournal(link, "worker-a")
        @test islink(link)
        root = joinpath(directory, "journal")
        journal = ResourceJournal(root, "worker-a")
        receipt = reserve_resource!(journal, journal_fence(), :podman, repeat("e",64))
        receipt_path = joinpath(root, string(receipt.id) * ".json")
        original = read(receipt_path, String)
        close(journal)
        pending = joinpath(root, ".pending-" * string(uuid4()))
        write(pending, "partial")
        chmod(pending, 0o600)
        write(receipt_path, "{corrupt")
        @test_throws ArgumentError ResourceJournal(root, "worker-a")
        @test read(receipt_path, String) == "{corrupt"
        @test isfile(pending) # Validate all committed files before removing fragments.
        rm(receipt_path)
        @test ccall(:mkfifo, Cint, (Cstring, Cuint), receipt_path, 0o600) == 0
        @test_throws ArgumentError ResourceJournal(root, "worker-a")
        @test ispath(pending)
        rm(receipt_path)
        write(receipt_path, original)
        chmod(receipt_path, 0o644)
        @test_throws ArgumentError ResourceJournal(root, "worker-a")
        @test stat(receipt_path).mode & 0o777 == 0o644
        chmod(receipt_path, 0o600)
        hardlink = joinpath(directory, "same-inode")
        Base.Filesystem.hardlink(receipt_path, hardlink)
        @test_throws ArgumentError ResourceJournal(root, "worker-a")
        @test isfile(hardlink)
        rm(hardlink)
        unknown = joinpath(root, "foreign-resource")
        write(unknown, "do not delete")
        @test_throws ArgumentError ResourceJournal(root, "worker-a")
        @test read(unknown, String) == "do not delete"
        @test isfile(pending)
        rm(unknown)
        restored = ResourceJournal(root, "worker-a")
        try
            @test !ispath(pending)
            @test length(resource_receipts(restored)) == 1
            # Replacing a held lock cannot allow the old owner to keep writing.
            oldlock = joinpath(directory, "held-lock")
            Base.Filesystem.rename(joinpath(root, "owner.lock"), oldlock)
            write(joinpath(root, "owner.lock"), "")
            chmod(joinpath(root, "owner.lock"), 0o600)
            @test_throws ArgumentError reserve_resource!(restored, journal_fence(), :podman, repeat("e",64))
        finally
            close(restored)
        end
    end
end

@testset "kernel ownership is released on agent death, receipts are not leases" begin
    mktempdir() do directory
        root = joinpath(directory, "journal")
        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        child = joinpath(@__DIR__, "journal_child.jl")
        project = normpath(joinpath(@__DIR__, ".."))
        process = open(`$julia --startup-file=no --compiled-modules=existing --project=$project $child $root`, "r+")
        reader = @async readline(process)
        try
            @test timedwait(() -> istaskdone(reader), 15; pollint=0.01) == :ok
            @test fetch(reader) == "journal-ready"
            @test_throws ArgumentError ResourceJournal(root, "worker-a")
            kill(process, Base.SIGKILL)
            @test timedwait(() -> !Base.process_running(process), 5; pollint=0.01) == :ok
            wait(process, false)
            restored = ResourceJournal(root, "worker-a")
            try
                receipt = only(resource_receipts(restored))
                @test receipt.physical_id === nothing
                @test receipt.fence.worker_id == "worker-a"
                @test length(readdir(root)) == 3
            finally
                close(restored)
            end
        finally
            Base.process_running(process) && kill(process, Base.SIGKILL)
            wait(process, false)
            close(process)
            @test timedwait(() -> istaskdone(reader), 5; pollint=0.01) == :ok
        end
    end
end
