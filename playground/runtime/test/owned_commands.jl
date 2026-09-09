using Test
using LineCableModelsRuntime
const CommandRT = LineCableModelsRuntime

function fixture_command(mode)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    return setenv(`$julia --startup-file=no --history-file=no $(joinpath(@__DIR__, "command_child.jl")) $mode`,
        Dict("PATH"=>get(ENV, "PATH", "/usr/bin:/bin")))
end
command_error(action) = try action(); nothing catch error; error end

@testset "host command lifetime, output and authority bounds" begin
    @test_throws ArgumentError CommandRunner(capacity=0)
    @test_throws ArgumentError CommandRunner(capacity=true)
    @test_throws ArgumentError CommandRunner(timeout_seconds=Inf)
    @test_throws ArgumentError CommandRunner(maximum_bytes=0)
    runner = CommandRunner(maximum_bytes=128, cleanup_seconds=0.2)
    try
        @test isempty(runner.active)
        result = run_owned_command!(runner, fixture_command("echo"))
        @test result.exitcode == 17
        @test result.output == "private-stdout"
        @test result.diagnostic == "private-stderr"
        @test !occursin("private-stdout", repr(MIME"text/plain"(), result))
        @test isempty(runner.active)
        result = run_owned_command!(runner, fixture_command("stdin"))
        @test result.output == "closed"
        @test isempty(runner.active)
        task = @async run_owned_command!(runner, fixture_command("flood"))
        @test timedwait(() -> istaskdone(task), 10; pollint=0.01) == :ok
        error = command_error(() -> fetch(task))
        @test error isa TaskFailedException
        @test occursin("output_limit", sprint(showerror, error))
        @test !occursin("private-flood", sprint(showerror, error))
        @test isempty(runner.active)
        error = command_error(() -> run_owned_command!(runner, fixture_command("wait"); timeout_seconds=0.5))
        @test error isa CommandFailure && error.code == :deadline
        @test isempty(runner.active)
        error = command_error(() -> run_owned_command!(runner, `/missing/private-command secret-argument`))
        @test error isa CommandFailure && error.code == :command_failed
        @test !occursin("secret-argument", sprint(showerror, error))
        @test isempty(runner.active)
        token = CommandRT.ExecutionCore.CancellationToken()
        CommandRT.ExecutionCore.cancel!(token)
        error = command_error(() -> run_owned_command!(runner, fixture_command("wait"); token))
        @test error isa CommandFailure && error.code == :canceled
        @test isempty(runner.active)
    finally
        close(runner)
    end
    @test close(runner) === nothing
    error = command_error(() -> run_owned_command!(runner, fixture_command("echo")))
    @test error isa CommandFailure && error.code == :closed
end

@testset "independent commands, bounded admission and exact cancellation" begin
    runner = CommandRunner(capacity=2, cleanup_seconds=0.2)
    token = CommandRT.ExecutionCore.CancellationToken()
    long = @async run_owned_command!(runner, fixture_command("wait"); token)
    try
        @test timedwait(() -> length(runner.active) == 1, 5; pollint=0.01) == :ok
        @test run_owned_command!(runner, fixture_command("echo")).exitcode == 17
        other = @async run_owned_command!(runner, fixture_command("wait"))
        @test timedwait(() -> length(runner.active) == 2, 5; pollint=0.01) == :ok
        error = command_error(() -> run_owned_command!(runner, fixture_command("echo")))
        @test error isa CommandFailure && error.code == :busy
        CommandRT.ExecutionCore.cancel!(token)
        @test timedwait(() -> istaskdone(long), 5; pollint=0.01) == :ok
        @test occursin("canceled", sprint(showerror, command_error(() -> fetch(long))))
        @test !istaskdone(other)
        @test length(runner.active) == 1
        @test close(runner) === nothing
        @test timedwait(() -> istaskdone(other), 5; pollint=0.01) == :ok
        @test isempty(runner.active)
    finally
        close(runner)
    end
end

@testset "host command unresolved cleanup remains owned and retryable" begin
    runner = CommandRunner(cleanup_seconds=0.05)
    release = Channel{Nothing}(1)
    task = @async run_owned_command!(runner, fixture_command("wait"))
    try
        @test timedwait(() -> length(runner.active) == 1, 5; pollint=0.01) == :ok
        handle = only(values(runner.active))
        # Fault injection: an independently blocked owned reader cannot be forgotten.
        push!(handle.readers, @async take!(release))
        error = command_error(() -> close(runner))
        @test error isa CommandFailure && error.code == :cleanup_unresolved
        @test runner.closed
        @test get(runner.active, handle.id, nothing) === handle
        @test !Base.process_running(handle.process)
        @test timedwait(() -> istaskdone(task), 5; pollint=0.01) == :ok
        @test occursin("cleanup_unresolved", sprint(showerror, command_error(() -> fetch(task))))
        put!(release, nothing)
        @test close(runner) === nothing
        @test isempty(runner.active)
    finally
        isready(release) || put!(release, nothing)
        close(runner)
    end
end

@testset "container command environment does not forward private authority" begin
    withenv("NATS_CONNECT_URL"=>"secret-broker", "AWS_SECRET_ACCESS_KEY"=>"secret-storage",
            "LCM_PROXY_KEY"=>"secret-proxy", "DOCKER_HOST"=>"tcp://remote:2375",
            "DOCKER_CONTEXT"=>"remote", "CONTAINER_HOST"=>"ssh://remote/run/podman.sock",
            "CONTAINER_CONNECTION"=>"remote") do
        env = CommandRT.container_command_environment()
        @test env["LC_ALL"] == "C"
        for key in ("NATS_CONNECT_URL", "AWS_SECRET_ACCESS_KEY", "LCM_PROXY_KEY",
                "DOCKER_HOST", "DOCKER_CONTEXT", "CONTAINER_HOST", "CONTAINER_CONNECTION")
            @test !haskey(env, key)
        end
    end
end
