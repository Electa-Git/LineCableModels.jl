using Test
using LineCableModelsRuntime
const TerminalRT=LineCableModelsRuntime
terminal_error(f)=try f();nothing catch error;error end

function terminal_fixture(mode;limits=TerminalRT.TerminalIOLimits(lifetime_seconds=45))
    process=TerminalRT.TerminalProcess(;limits)
    julia=joinpath(Sys.BINDIR,Base.julia_exename())
    command=setenv(`$julia --startup-file=no --history-file=no --color=yes $(joinpath(@__DIR__,"terminal_child.jl")) $mode`,
        Dict("TERM"=>"xterm-256color","JULIA_NUM_THREADS"=>"1","OPENBLAS_NUM_THREADS"=>"1"))
    try
        TerminalRT.start_terminal!(process,command)
    catch
        close(process);rethrow()
    end
    return process
end

function terminal_text(process)
    return lock(process.lock) do
        String(TerminalRT.read_terminal(process.output,0,length(process.output.data)).bytes)
    end
end
terminal_write(process,text)=TerminalRT.write_terminal!(process,Vector{UInt8}(codeunits(text)))
function terminal_expect(process,text;seconds=15)
    result=timedwait(()->occursin(text,terminal_text(process)) || process.cleanup_complete,seconds;pollint=0.01)
    @test result==:ok
    matched=occursin(text,terminal_text(process))
    @test matched
end

@testset "terminal byte bounds, cursors and private display" begin
    for options in ((chunk_bytes=true,), (chunk_bytes=16385,), (input_bytes=1,),
            (output_bytes=1,), (output_bytes_per_second=1,), (write_seconds=Inf,),
            (lifetime_seconds=0,), (cleanup_seconds=true,))
        @test_throws ArgumentError TerminalRT.TerminalIOLimits(;options...)
    end
    buffer=TerminalRT.TerminalBuffer(7)
    TerminalRT.append_terminal!(buffer,UInt8[1,2,3,4])
    first=TerminalRT.read_terminal(buffer,0,3)
    @test first.bytes==UInt8[1,2,3] && first.cursor==3 && first.sequence==4 && !first.gap
    TerminalRT.append_terminal!(buffer,UInt8[5,6,7,8,9,10])
    lagged=TerminalRT.read_terminal(buffer,0,4)
    @test lagged.gap && lagged.bytes==UInt8[4,5,6,7] && lagged.cursor==7
    @test TerminalRT.read_terminal(buffer,7,4).bytes==UInt8[8,9,10]
    @test isempty(TerminalRT.read_terminal(buffer,10,4).bytes)
    @test_throws ArgumentError TerminalRT.read_terminal(buffer,-1,1)
    @test_throws ArgumentError TerminalRT.read_terminal(buffer,11,1)
    @test_throws ArgumentError TerminalRT.read_terminal(buffer,true,1)
    @test_throws ArgumentError TerminalRT.read_terminal(buffer,0,0)
    TerminalRT.append_terminal!(buffer,Vector{UInt8}(codeunits("private-terminal-output")))
    @test length(buffer.data)==7
    @test !occursin("output",repr(buffer))
    @test !occursin("output",repr(TerminalRT.read_terminal(buffer,0,7)))
    buffer.sequence=9007199254740991
    @test_throws TerminalRT.TerminalFailure TerminalRT.append_terminal!(buffer,UInt8[1])
    process=TerminalRT.TerminalProcess()
    @test process.process===nothing && process.master == -1 && process.pump===nothing
    @test_throws ArgumentError TerminalRT.start_terminal!(process,`true`)
    @test_throws ArgumentError TerminalRT.resize_terminal!(process,0,20)
    @test_throws TerminalRT.TerminalFailure terminal_write(process,"private-input")
    @test close(process)===nothing
    @test process.cleanup_complete
    @test_throws TerminalRT.TerminalFailure TerminalRT.start_terminal!(process,setenv(`true`,Dict{String,String}()))
end

@testset "real PTY transport ownership and scheduler independence" begin
    process=terminal_fixture("echo")
    try
        terminal_expect(process,"READY")
        @test ccall(:isatty,Cint,(Cint,),process.master)==1
        @test ccall(:fcntl,Cint,(Cint,Cint),process.master,3) & Base.JL_O_NONBLOCK != 0
        @test ccall(:fcntl,Cint,(Cint,Cint),process.master,1) & 1 == 1
        terminal_write(process,"private-input-λ\n")
        terminal_expect(process,"ECHO:private-input-λ")
        @test !occursin("private-input",repr(MIME"text/plain"(),process))
        @test_throws ArgumentError terminal_write(process,repeat("x",process.limits.chunk_bytes+1))
        @test_throws TerminalRT.TerminalFailure TerminalRT.start_terminal!(process,setenv(`true`,Dict{String,String}()))
        terminal_write(process,"quit\n")
        @test timedwait(()->process.cleanup_complete,10;pollint=0.01)==:ok
        @test process.process.exitcode==0 && process.master == -1
        @test istaskdone(process.pump)
    finally
        close(process)
    end
    @test close(process)===nothing
    failed=TerminalRT.TerminalProcess()
    error=terminal_error(()->TerminalRT.start_terminal!(failed,setenv(`/missing/private-command secret-input`,Dict{String,String}())))
    @test error isa TerminalRT.TerminalFailure && error.code==:spawn_failed
    @test !occursin("secret-input",sprint(showerror,error))
    @test failed.cleanup_complete && failed.master == -1
    close(failed)

    flood=terminal_fixture("flood";limits=TerminalRT.TerminalIOLimits(output_bytes_per_second=8192,lifetime_seconds=10))
    other=terminal_fixture("echo")
    try
        terminal_expect(other,"READY")
        terminal_write(other,"survives\n");terminal_expect(other,"ECHO:survives")
        @test timedwait(()->flood.cleanup_complete,10;pollint=0.01)==:ok
        @test flood.failure==:output_rate_limit
        @test length(flood.output.data)==flood.limits.output_bytes
        @test !process_running(flood.process) && flood.master == -1
        @test !other.closing
    finally
        close(flood);close(other)
    end

    stalled=terminal_fixture("stall";limits=TerminalRT.TerminalIOLimits(write_seconds=0.2,lifetime_seconds=10))
    try
        terminal_expect(stalled,"READY")
        block=fill(UInt8('x'),stalled.limits.chunk_bytes)
        lock(stalled.lock) do
            for _ in 1:div(stalled.limits.input_bytes,length(block))
                TerminalRT.write_terminal!(stalled,block)
            end
            @test_throws TerminalRT.TerminalFailure TerminalRT.write_terminal!(stalled,block)
            @test stalled.pending_bytes==stalled.limits.input_bytes
        end
        @test timedwait(()->stalled.cleanup_complete,10;pollint=0.01)==:ok
        @test stalled.failure==:input_stalled
        @test stalled.pending_bytes==0 && isempty(stalled.input)
    finally
        close(stalled)
    end
    expires=terminal_fixture("stall";limits=TerminalRT.TerminalIOLimits(lifetime_seconds=0.5))
    try
        @test timedwait(()->expires.cleanup_complete,10;pollint=0.01)==:ok
        @test expires.failure==:lifetime_limit
    finally
        close(expires)
    end
end

@testset "terminal partial cleanup retains ownership for an exact retry" begin
    process=terminal_fixture("echo";limits=TerminalRT.TerminalIOLimits(cleanup_seconds=0.05,lifetime_seconds=10))
    gate=Channel{Nothing}(1)
    original=process.pump
    try
        terminal_expect(process,"READY")
        # Inject a stuck owned task, retaining the actual I/O pump separately.
        process.pump=@async take!(gate)
        fd=process.master
        error=terminal_error(()->close(process))
        @test error isa TerminalRT.TerminalFailure && error.code==:cleanup_unresolved
        @test process.master==fd && !process.cleanup_complete
        @test !process_running(process.process)
        @test timedwait(()->istaskdone(original),2;pollint=0.01)==:ok
        put!(gate,nothing)
        @test close(process)===nothing
        @test process.cleanup_complete && process.master == -1 && istaskdone(process.pump)
        @test_throws TerminalRT.TerminalFailure TerminalRT.resize_terminal!(process,80,24)
    finally
        isready(gate) || istaskdone(process.pump) || put!(gate,nothing)
        close(process)
    end
    # A closed descriptor must not accumulate across repeated short-lived PTYs.
    baseline=length(readdir("/proc/self/fd"))
    for _ in 1:8
        child=TerminalRT.TerminalProcess()
        try
            TerminalRT.start_terminal!(child,setenv(`/usr/bin/true`,Dict{String,String}()))
            @test timedwait(()->child.cleanup_complete,5;pollint=0.01)==:ok
            @test child.failure===nothing
        finally
            close(child)
        end
    end
    @test length(readdir("/proc/self/fd"))<=baseline
end

@testset "actual Julia REPL, not a substitute line parser" begin
    process=terminal_fixture("repl")
    try
        terminal_expect(process,"julia>")
        terminal_write(process,"α = 41\rprintln(\"RESULT=\",α+1)\r")
        terminal_expect(process,"RESULT=42")
        private=terminal_fixture("repl")
        try
            terminal_expect(private,"julia>")
            terminal_write(private,"println(\"PRIVATE=\",isdefined(Main,:α))\r")
            terminal_expect(private,"PRIVATE=false")
            @test private.id != process.id
        finally
            close(private)
        end
        terminal_write(process,"function fixture_double(x)\r2x\rend\rprintln(\"MULTILINE=\",fixture_double(21))\r")
        terminal_expect(process,"MULTILINE=42")
        terminal_write(process,"println(\"COMPLETE=\",fixture_dou\t(21))\r")
        terminal_expect(process,"COMPLETE=42")
        TerminalRT.resize_terminal!(process,120,40)
        terminal_write(process,"println(\"SIZE=\",displaysize(stdout))\r")
        terminal_expect(process,"SIZE=(40, 120)")
        terminal_write(process,"println(\"HISTORY=\",42)\r")
        terminal_expect(process,"HISTORY=42")
        terminal_write(process,"\e[A\r")
        @test timedwait(()->length(findall("HISTORY=42",terminal_text(process)))>=2,10;pollint=0.01)==:ok
        terminal_write(process,"while true; sleep(0.1); end\r")
        sleep(0.2)
        terminal_write(process,"\x03")
        terminal_expect(process,"InterruptException")
        terminal_write(process,"println(\"AFTER-INTERRUPT=\",α+1)\r")
        terminal_expect(process,"AFTER-INTERRUPT=42")
        terminal_write(process,"exit()\r")
        @test timedwait(()->process.cleanup_complete,10;pollint=0.01)==:ok
        @test process.process.exitcode==0
    finally
        close(process)
    end
end
