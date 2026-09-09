function blocked_executor(mode="blocked"; redirected=false)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    child = joinpath(@__DIR__, "blocked_executor.jl")
    command = setenv(`setpriv --pdeathsig KILL $julia --startup-file=no --history-file=no $child $mode`,
        ["PATH"=>ENV["PATH"], "JULIA_NUM_THREADS"=>"1", "OPENBLAS_NUM_THREADS"=>"1"])
    return ExecutorSupervisor(;command,discard_stderr=redirected)
end

@testset "redirected attached command retains original process ownership" begin
    @test_throws TypeError ExecutorSupervisor(;command=pipeline(`true`,`true`))
    supervisor = blocked_executor(;redirected=true)
    try
        process = start_executor!(supervisor,context())
        @test process isa Base.Process && process_running(process)
        @test supervisor.process === process && supervisor.generation == 1
        @test stop_executor!(supervisor;grace_seconds=0.1,kill_seconds=1) === nothing
        @test !process_running(process) && supervisor.process === nothing
        @test isempty(supervisor.io_tasks)
    finally
        stop_executor!(supervisor)
    end
end

@testset "bounded executor pipe submission and joined cleanup" begin
    for invalid in (0, -1, Inf, true, 11)
        @test_throws ArgumentError stop_executor!(ExecutorSupervisor(; command=`false`); grace_seconds=invalid)
        @test_throws ArgumentError stop_executor!(ExecutorSupervisor(; command=`false`); kill_seconds=invalid)
    end
    idle = ExecutorSupervisor(; command=`false`)
    canceled = context()
    cancel!(canceled.token)
    @test_throws E.OperationCanceled start_executor!(idle, canceled)
    @test idle.process === nothing && idle.generation == 0
    inputs = Dict{String,Any}("fill"=>repeat("x", 512 * 1024))
    spec = OperationSpec("fixture.blocked", identity, (_, p) -> p; timeout_seconds=0.15, execution_mode=:supervised)
    supervisor = blocked_executor()
    ticks = Ref(0)
    monitoring = Ref(true)
    monitor = @async while monitoring[]
        ticks[] += 1
        sleep(0.01)
    end
    try
        process = start_executor!(supervisor, context())
        elapsed = @elapsed @test_throws PermanentOperationError execute_supervised!(supervisor, spec, context(), inputs)
        @test elapsed < 9
        @test ticks[] > 5
        @test !process_running(process)
        @test supervisor.process === nothing
        @test isempty(supervisor.io_tasks)
    finally
        monitoring[] = false
        wait(monitor)
        stop_executor!(supervisor)
    end
    supervisor = blocked_executor()
    timer = nothing
    try
        process = start_executor!(supervisor, context())
        canceled = context()
        timer = Timer(0.15) do _
            cancel!(canceled.token)
        end
        patient = OperationSpec("fixture.blocked", identity, (_, p) -> p; timeout_seconds=60, execution_mode=:supervised)
        @test_throws E.OperationCanceled execute_supervised!(supervisor, patient, canceled, inputs)
        @test !process_running(process) && supervisor.process === nothing
        @test isempty(supervisor.io_tasks)
    finally
        timer === nothing || close(timer)
        stop_executor!(supervisor)
    end
    supervisor = blocked_executor()
    try
        process = start_executor!(supervisor, context())
        # Julia's signal-wait thread can handle TERM even when libc signal()
        # requests SIG_IGN. Stop this exact owned child instead; TERM remains
        # pending, so only the supervisor's KILL escalation can retire it.
        # Base exposes TERM/KILL but not STOP; obtain the platform's number.
        stop_signal = parse(Int, strip(read(`bash -c "kill -l STOP"`, String)))
        kill(process, stop_signal)
        stopped() = any(line -> startswith(line, "State:") && occursin("T (stopped)", line),
            readlines("/proc/$(getpid(process))/status"))
        @test timedwait(stopped, 2; pollint=0.025) == :ok
        elapsed = @elapsed @test stop_executor!(supervisor; grace_seconds=0.1, kill_seconds=1) === nothing
        @test elapsed < 3
        @test !process_running(process) && process.termsignal == Base.SIGKILL
        @test isempty(supervisor.io_tasks)
        @test stop_executor!(supervisor) === nothing
    finally
        stop_executor!(supervisor)
    end
    supervisor = blocked_executor("malformed")
    try
        process = start_executor!(supervisor, context())
        malformed = OperationSpec("fixture.malformed", identity, (_, p) -> p; timeout_seconds=10, execution_mode=:supervised)
        @test_throws RetryableOperationError execute_supervised!(supervisor, malformed, context(), Dict{String,Any}())
        @test !process_running(process) && supervisor.process === nothing
        @test isempty(supervisor.io_tasks)
    finally
        stop_executor!(supervisor)
    end
end

@testset "unresolved executor cleanup retains its original handle" begin
    supervisor = blocked_executor()
    gate = Channel{Nothing}(1)
    pending = nothing
    try
        process = start_executor!(supervisor, context())
        # Model a pipe task whose close has not completed; retain it as owned
        # until its completion is observed instead of releasing the slot early.
        pending = E.executor_io_task(() -> take!(gate), supervisor, process)
        @test_throws RetryableOperationError stop_executor!(supervisor; grace_seconds=0.1, kill_seconds=0.1)
        @test supervisor.process === process
        @test !process_running(process)
        @test pending in supervisor.io_tasks && !istaskdone(pending)
        put!(gate, nothing)
        wait(pending)
        first = @async stop_executor!(supervisor; grace_seconds=0.1, kill_seconds=0.1)
        second = @async stop_executor!(supervisor; grace_seconds=0.1, kill_seconds=0.1)
        @test fetch(first) === nothing && fetch(second) === nothing
        @test supervisor.process === nothing && isempty(supervisor.io_tasks)
    finally
        pending === nothing || istaskdone(pending) || (put!(gate, nothing); wait(pending))
        stop_executor!(supervisor)
    end
end
