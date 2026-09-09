using Test, UUIDs, LineCableModelsRuntime
const NativeSystemdRT = LineCableModelsRuntime

function test_native_recovery_systemd()
    directory = mktempdir(;prefix="lcm-native-recovery-",cleanup=false)
    worker = "audit-" * string(uuid4())
    journal = ResourceJournal(joinpath(directory,"journal"),worker;capacity=3)
    runner = CommandRunner()
    passed = false
    try
        scope = native_scope(runner)
        @test occursin(r"^[a-f0-9]{64}$",scope)
        receipts = ResourceReceipt[]
        function reserve()
            fence = NativeSystemdRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"audit","cleanup",
                worker,string(uuid4()),string(uuid4()),"cleanup-fixture","1.0.0",repeat("a",64),1)
            reserve_resource!(journal,fence,:native,scope)
        end
        # Warm the absent-resource path before any real service exists.
        absent = reserve()
        @test remove_owned_native!(journal,runner,absent;timeout_seconds=30)
        function launch_fixture(receipt;properties=String[],nonblocking=false)
            # Harmless, finite cleanup fixture, not the scientific launch path:
            # no numerical package, terminal, CPU-admission bypass or user code.
            command = [Sys.which("systemd-run"),"--user","--quiet","--unit=" * resource_name(receipt),
                "--slice=app.slice","--property=Description=" * native_resource_description(receipt),
                "--property=Type=exec","--property=ExitType=main","--property=RemainAfterExit=no",
                "--property=KillMode=control-group","--property=KillSignal=SIGINT","--property=SendSIGKILL=yes",
                "--property=FinalKillSignal=SIGKILL","--property=TimeoutStopSec=2","--property=TimeoutStopFailureMode=terminate",
                "--property=Restart=no","--property=RuntimeMaxSec=120","--property=MemoryMax=67108864",
                "--property=TasksMax=8","--property=StandardInput=null","--property=StandardOutput=null",
                "--property=StandardError=null"]
            nonblocking && push!(command,"--no-block")
            append!(command,["--property=" * property for property in properties])
            append!(command,["--",Sys.which("sleep"),"120"])
            result = run_owned_command!(runner,setenv(Cmd(command),NativeSystemdRT.container_command_environment()))
            if result.exitcode != 0
                diagnostic = joinpath(directory,"launch-error.txt")
                open(diagnostic,"w") do io; chmod(diagnostic,0o600); write(io,result.diagnostic); end
                error("Owned native cleanup fixture launch failed.")
            end
        end
        for _ in 1:2
            receipt = reserve()
            push!(receipts,receipt)
            launch_fixture(receipt)
        end
        first,second = receipts
        one,two = inspect_native_unit(runner,first),inspect_native_unit(runner,second)
        @test one !== nothing && two !== nothing
        @test one.pid > 0 && two.pid > 0 && one.pid != two.pid
        @test strip(read("/proc/$(one.pid)/cgroup",String)) == "0::" * one.cgroup
        @test !NativeSystemdRT.native_group_empty(one.cgroup)
        @test all(r->r.physical_id===nothing,resource_receipts(journal))
        @test remove_owned_native!(journal,runner,first;timeout_seconds=30)
        @test inspect_native_unit(runner,first) === nothing
        @test NativeSystemdRT.native_group_empty(one.cgroup)
        survivor = inspect_native_unit(runner,second)
        @test survivor.invocation == two.invocation && survivor.pid == two.pid
        @test only(resource_receipts(journal)).id == second.id
        # Simulate a parent restart: close/reopen the persistent ownership file
        # while the second exact service remains live under systemd. Use the
        # actual root post-stop dispatcher, not just the native helper.
        close(journal)
        try
            @test recover_agent_resources!(joinpath(directory,"journal"),worker) === nothing
        finally
            journal = ResourceJournal(joinpath(directory,"journal"),worker;capacity=3)
        end
        @test inspect_native_unit(runner,second) === nothing
        @test NativeSystemdRT.native_group_empty(two.cgroup)
        @test isempty(resource_receipts(journal))
        @test recover_native!(journal,runner) === nothing
        # Hold one owned start job briefly. Its dependent has not entered
        # service_start yet, so it has no invocation or process to bind.
        barrier,queued = reserve(),reserve()
        launch_fixture(barrier;properties=["ExecStartPost=" * Sys.which("sleep") * " 30"],nonblocking=true)
        launch_fixture(queued;properties=["After=" * resource_name(barrier)],nonblocking=true)
        waiting = inspect_native_unit(runner,queued)
        blocker = inspect_native_unit(runner,barrier)
        @test waiting !== nothing && waiting.invocation === nothing && waiting.pid == 0
        @test !NativeSystemdRT.native_jobs_clear(runner,queued)
        @test remove_owned_native!(journal,runner,queued;timeout_seconds=30)
        @test inspect_native_unit(runner,queued) === nothing && NativeSystemdRT.native_group_empty(waiting.cgroup)
        @test inspect_native_unit(runner,barrier).invocation == blocker.invocation
        @test remove_owned_native!(journal,runner,barrier;timeout_seconds=30)
        @test isempty(resource_receipts(journal)) && NativeSystemdRT.native_group_empty(blocker.cgroup)
        passed = true
    finally
        try
            recover_native!(journal,runner;timeout_seconds=30)
        finally
            close(journal);close(runner)
            println("Native recovery diagnostics: ",directory)
        end
    end
    passed || error("Native service recovery test did not complete.")
    println("All exact harmless native services and the never-started queued unit removed; no executor admission claimed.")
end

@testset "actual native transient service cleanup and restart recovery" begin
    test_native_recovery_systemd()
end
