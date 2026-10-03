@testmodule PSCADReaderFailures begin
    const RemoteConfig = NamedTuple
    const fail_at_eof = Ref(false)
    const readers_started = Channel{Nothing}(2)
    const child = Ref{Base.Process}()

    remote_command(config::RemoteConfig, ::AbstractString) = config.command
    function run(command; wait)
        child[] = Base.run(command; wait)
    end

    struct Lines{T}
        source::T
    end
    function eachline(pipe::Pipe)
        put!(readers_started, nothing)
        return Lines(Base.eachline(pipe))
    end
    function Base.iterate(lines::Lines, state...)
        entry = iterate(lines.source, state...)
        entry === nothing && fail_at_eof[] && error("synthetic reader cleanup failure")
        entry === nothing || first(entry) != "FAIL" || error("synthetic reader failure")
        return entry
    end

    # Execute the production transport method with isolated IO fault injection.
    Base.include(@__MODULE__,
        joinpath(@__DIR__, "../../../ext/LineCableModelsPSCADExt/remote/remote.jl")) do expression
        expression isa Expr && expression.head === :function &&
        expression.args[1] isa Expr && expression.args[1].head === :call &&
        expression.args[1].args[1] === :_run_remote ? expression : nothing
    end
end

@testitem "PSCAD / reader failures propagate and cancellation retains secondary diagnostics" tags=[:integration, :pscad] default_imports=false setup=[PSCADReaderFailures] begin
    using Test, Logging
    const P=PSCADReaderFailures
    for stream in ("stdout", "stderr")
        script="println($stream, \"partial\"); println($stream, \"FAIL\"); flush($stream); sleep(30)"
        command=`$(Base.julia_cmd()) --startup-file=no --handle-signals=no --project=@stdlib -e $script`
        mktempdir() do root
            stdout_path, stderr_path=joinpath(root, "out.log"), joinpath(root, "err.log")
            caught=try
                P._run_remote((; command, timeout = 10.0), "unused"; stdout_path, stderr_path)
            catch exception
                exception
            end
            @test caught isa TaskFailedException
            @test occursin("synthetic reader failure", sprint(showerror, caught))
            @test process_exited(P.child[])
            @test read(stream == "stdout" ? stdout_path : stderr_path, String) ==
                  "partial\n"
            take!(P.readers_started)
            take!(P.readers_started)
        end
    end

    command=`$(Base.julia_cmd()) --startup-file=no --handle-signals=no --project=@stdlib -e 'sleep(30)'`
    P.fail_at_eof[]=true
    logger=Test.TestLogger(min_level = Logging.Warn)
    canceled=Ref(false)
    interrupted=InterruptException()
    try
        task=with_logger(logger) do
            @async try
                P._run_remote((; command, timeout = 10.0), "unused";
                    on_interrupt = () -> (canceled[] = true))
            catch exception
                exception
            end
        end
        take!(P.readers_started)
        take!(P.readers_started)
        schedule(task, interrupted; error = true)
        @test fetch(task) === interrupted
        @test canceled[]
        @test process_exited(P.child[])
        @test Set(record.kwargs[:stream] for record in logger.logs) ==
              Set((:stdout, :stderr))
        @test all(record -> occursin("reader failed during cleanup", string(record.message)), logger.logs)
        @test all(
            record -> occursin("synthetic reader cleanup failure",
                sprint(showerror, first(record.kwargs[:exception]))),
            logger.logs)
    finally
        P.fail_at_eof[]=false
    end
end
