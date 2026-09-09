using Test, UUIDs
using LineCableModelsRuntime
const StoppedRT = LineCableModelsRuntime

function test_stopped_container_recovery()
image = get(ENV, "LCM_TEST_STOPPED_CONTAINER_IMAGE", "")
occursin(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*@sha256:[a-f0-9]{64}$", image) ||
    error("Set LCM_TEST_STOPPED_CONTAINER_IMAGE to an already cached digest-pinned image.")
# This gate deliberately tests acquisition/recovery only. No start/run command is
# issued; every configured limit remains present even when this host cannot run it.
runner = CommandRunner()
test_directory = mktempdir(; prefix="lcm-stopped-recovery-", cleanup=false)
journal = nothing
passed = false
cleaned = false
try
    host = check_container_host(runner; requested=get(ENV, "LCM_TEST_CONTAINER_RUNTIME", "auto"))
    inspect_image = StoppedRT.scoped_container_command(runner, host,
        ["image", "inspect", "--format", "{{json .RepoDigests}}", image])
    inspect_image.exitcode == 0 || error("The test image must already exist locally; no pull is permitted.")
    image in StoppedRT.JSON3.read(inspect_image.output) || error("The cached image digest does not match.")
    scope = container_scope(runner, host)
    before = StoppedRT.scoped_container_command(runner, host,
        ["container", "ls", "--all", "--no-trunc", "--format", "{{.ID}}"])
    before.exitcode == 0 || error("Initial container inventory is unavailable.")
    journal = ResourceJournal(joinpath(test_directory, "journal"), "stopped-audit"; capacity=1)
    fence = StoppedRT.Protocol.AssignmentFence(string(uuid4()), string(uuid4()), "audit", "cleanup",
        "stopped-audit", string(uuid4()), string(uuid4()), "cleanup-audit", "1.0.0",
        last(split(image, "@sha256:")), 1)
    receipt = reserve_resource!(journal, fence, host.engine.name, scope)
    arguments = ["container", "create", "--name", resource_name(receipt),
        "--pull=never", "--read-only", "--network=none", "--user=1000:1000",
        "--cap-drop=ALL", "--security-opt=no-new-privileges", "--log-driver=none",
        "--restart=no", "--cpus=1", "--memory=67108864", "--memory-swap=67108864",
        "--pids-limit=16", "--shm-size=65536",
        "--tmpfs", "/tmp:rw,nosuid,nodev,noexec,size=8388608,mode=1777",
        "--entrypoint=/usr/bin/true"]
    host.engine.name == :podman && append!(arguments, ["--read-only-tmpfs=false", "--image-volume=ignore"])
    for (key, value) in sort!(collect(resource_labels(receipt)); by=first)
        append!(arguments, ["--label", key * "=" * value])
    end
    push!(arguments, image)
    created = StoppedRT.scoped_container_command(runner, host, arguments)
    if created.exitcode != 0
        diagnostic = joinpath(test_directory, "create-diagnostic.txt")
        write(diagnostic, created.diagnostic)
        chmod(diagnostic, 0o600)
        error("Stopped acquisition failed; bounded private diagnostics were retained.")
    end
    id = strip(created.output)
    @testset "actual stopped-container acquisition and exact recovery" begin
        @test occursin(r"^[a-f0-9]{64}$", id)
        inspected = StoppedRT.scoped_container_command(runner, host,
            ["container", "inspect", "--format", "{{json .}}", id])
        @test inspected.exitcode == 0
        object = StoppedRT.JSON3.read(inspected.output)
        @test object.State.Running === false
        @test object.State.Pid == 0
        @test only(resource_receipts(journal)).physical_id === nothing
        # Leave the receipt at its pre-bind crash gap, then reopen its durable owner.
        close(journal)
        journal = ResourceJournal(joinpath(test_directory, "journal"), "stopped-audit"; capacity=1)
        recovered = only(resource_receipts(journal))
        @test recovered.id == receipt.id
        @test recovered.physical_id === nothing
        @test recover_containers!(journal, runner, host) === nothing
        @test isempty(resource_receipts(journal))
        @test remove_owned_container!(journal, runner, host, receipt)
        after = StoppedRT.scoped_container_command(runner, host,
            ["container", "ls", "--all", "--no-trunc", "--format", "{{.ID}}"])
        @test after.exitcode == 0
        @test sort(split(strip(after.output), '\n')) == sort(split(strip(before.output), '\n'))
        @test isempty(runner.active)
    end
    passed = true
    println("One newly created stopped container was removed by exact receipt-verified ID.")
    println("No container process was started; this is not an effective-limit or scientific-execution gate.")
finally
    try
        if journal !== nothing
            journal.closed && (journal = ResourceJournal(joinpath(test_directory, "journal"), "stopped-audit"; capacity=1))
            host_cleanup = check_container_host(runner; requested=get(ENV, "LCM_TEST_CONTAINER_RUNTIME", "auto"))
            recover_containers!(journal, runner, host_cleanup)
            cleaned = isempty(resource_receipts(journal))
        else
            cleaned = true
        end
    finally
        journal === nothing || close(journal)
        close(runner)
        if passed && cleaned
            realpath(test_directory) == test_directory &&
                startswith(basename(test_directory), "lcm-stopped-recovery-") ||
                error("Refusing unexpected test directory cleanup.")
            rm(test_directory; recursive=true)
        else
            println(stderr, "Owned recovery diagnostics retained: ", test_directory)
        end
    end
end
end

test_stopped_container_recovery()
