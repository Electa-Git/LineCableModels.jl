function agent_lease_fixture(; capacity=1, max_history=32)
    profiles = ProfileRegistry()
    register!(profiles, ProfileDefinition("line-parameters", "/approved/project", repeat("a", 64);
        operations=("system.echo",)))
    clock = Ref(0.0)
    ledger = AgentLeaseLedger("worker-a", profiles; capacity, max_history, clock=()->clock[])
    coordinator = string(uuid4())
    probe = RT.Protocol.WorkerProbe("2.0", "worker-a", coordinator, string(uuid4()))
    fence = RT.Protocol.AssignmentFence(string(uuid4()), string(uuid4()), "alice", "main", "worker-a",
        ledger.boot_id, coordinator, "line-parameters", "1.0.0", repeat("a", 64), 1)
    return (; ledger, clock, coordinator, probe, fence)
end
change_runtime_record(value::T; changes...) where T =
    T((get(changes, field, getfield(value, field)) for field in fieldnames(T))...)
lease_command(fence; action="grant", revision=1, duration_ms=5000, request_id=string(uuid4())) =
    RT.Protocol.LeaseControl("2.0", request_id, action, fence, revision, duration_ms)

@testset "agent authority uses its own clock and does not renew on duplicate delivery" begin
    (; ledger, clock, probe, fence) = agent_lease_fixture()
    grant = lease_command(fence)
    @test !handle_lease_control!(ledger, grant).accepted
    @test receive_probe!(ledger, probe)
    @test handle_lease_control!(ledger, grant).accepted
    @test agent_lease_usable(ledger, fence)
    @test !handle_lease_control!(ledger, lease_command(change_runtime_record(fence; worker_boot=string(uuid4())))).accepted
    original_expiry = only(values(ledger.leases)).expires_at
    clock[] = 2
    @test handle_lease_control!(ledger, grant).accepted
    @test only(values(ledger.leases)).expires_at == original_expiry
    @test !handle_lease_control!(ledger, change_runtime_record(grant; duration_ms=6000)).accepted
    renewal = lease_command(fence; action="renew", revision=2)
    @test handle_lease_control!(ledger, renewal).accepted
    @test only(values(ledger.leases)).expires_at == 7
    clock[] = 3
    @test handle_lease_control!(ledger, renewal).accepted
    @test only(values(ledger.leases)).expires_at == 7
    @test !handle_lease_control!(ledger, grant).accepted
    clock[] = 7
    @test expire_agent_leases!(ledger) == [fence]
    @test !agent_lease_usable(ledger, fence)
    other = change_runtime_record(fence; lease_id=string(uuid4()), run_id=string(uuid4()))
    @test !handle_lease_control!(ledger, lease_command(other)).accepted # closing still occupies capacity
    @test_throws AccessDenied complete_agent_cleanup!(ledger, other)
    @test complete_agent_cleanup!(ledger, fence) === nothing
    @test handle_lease_control!(ledger, lease_command(other)).accepted
    @test !agent_lease_usable(ledger, fence)
end

@testset "release acknowledges actual cleanup and fences delayed grants" begin
    (; ledger, probe, fence) = agent_lease_fixture()
    receive_probe!(ledger, probe)
    grant = lease_command(fence)
    handle_lease_control!(ledger, grant)
    @test_throws AccessDenied complete_agent_cleanup!(ledger, fence)
    release = lease_command(fence; action="release", revision=2, duration_ms=0)
    @test handle_lease_control!(ledger, release) === nothing
    @test !agent_lease_usable(ledger, fence)
    @test handle_lease_control!(ledger, release) === nothing
    ack = complete_agent_cleanup!(ledger, fence)
    @test ack.accepted && ack.request_id == release.request_id && ack.reason == "released"
    @test complete_agent_cleanup!(ledger, fence) == ack
    @test handle_lease_control!(ledger, release) == ack
    @test !handle_lease_control!(ledger, grant).accepted
    next_fence = change_runtime_record(fence; generation=2, lease_id=string(uuid4()))
    early_release = lease_command(next_fence; action="release", revision=2, duration_ms=0)
    @test handle_lease_control!(ledger, early_release).accepted
    @test !handle_lease_control!(ledger, lease_command(next_fence)).accepted
    newest = change_runtime_record(fence; generation=3, lease_id=string(uuid4()))
    @test handle_lease_control!(ledger, lease_command(newest)).accepted
    @test !handle_lease_control!(ledger, early_release).accepted
    @test agent_lease_usable(ledger, newest)
end

@testset "coordinator replacement must finish cleanup before advertising presence" begin
    (; ledger, clock, probe, fence) = agent_lease_fixture()
    receive_probe!(ledger, probe)
    handle_lease_control!(ledger, lease_command(fence; duration_ms=60_000))
    replacement = change_runtime_record(probe; coordinator_id=string(uuid4()), challenge=string(uuid4()))
    @test !receive_probe!(ledger, replacement)
    @test !agent_lease_usable(ledger, fence)
    @test_throws AccessDenied receive_probe!(ledger, probe) # old probe cannot cancel reconciliation
    @test expire_agent_leases!(ledger) == [fence]
    complete_agent_cleanup!(ledger, fence)
    @test receive_probe!(ledger, replacement)
    @test ledger.coordinator_id == replacement.coordinator_id
    @test_throws AccessDenied receive_probe!(ledger, probe)
    @test !handle_lease_control!(ledger, lease_command(fence)).accepted
    newer = change_runtime_record(fence; coordinator_id=replacement.coordinator_id,
        generation=2, lease_id=string(uuid4()))
    @test handle_lease_control!(ledger, lease_command(newer; duration_ms=60_000)).accepted
    clock[] = 11 # coordinator presence expires before this long grant
    @test !agent_lease_usable(ledger, newer)
    @test only(values(ledger.leases)).state == :closing
    @test receive_probe!(ledger, replacement) # restored heartbeat does not revive a revoked lease
    @test !agent_lease_usable(ledger, newer)
    complete_agent_cleanup!(ledger, newer)
    @test !agent_lease_usable(ledger, newer)
end

@testset "finite generation history fails closed instead of forgetting released authority" begin
    (; ledger, clock, probe, fence) = agent_lease_fixture(; max_history=1)
    receive_probe!(ledger, probe)
    release = lease_command(fence; action="release", revision=2, duration_ms=0)
    @test handle_lease_control!(ledger, release).accepted
    foreign = change_runtime_record(fence; run_id=string(uuid4()), lease_id=string(uuid4()))
    @test handle_lease_control!(ledger, lease_command(foreign)).reason == "history-capacity"
    newer = change_runtime_record(fence; generation=2, lease_id=string(uuid4()))
    @test handle_lease_control!(ledger, lease_command(newer)).accepted
    @test length(ledger.leases) == 1
    clock[] = -1
    @test !agent_lease_usable(ledger, newer) # even a broken/injected clock fails closed
    @test_throws ArgumentError AgentLeaseLedger("worker-a", ledger.profiles; capacity=true)
    @test_throws ArgumentError AgentLeaseLedger("worker-a", ledger.profiles; max_history=0)
    @test_throws ArgumentError AgentLeaseLedger("worker-a", ledger.profiles; presence_seconds=Inf)
end

