using Test

@testset "container entry fails before evaluation on missing/host isolation" begin
    project = normpath(joinpath(@__DIR__,".."))
    guard = normpath(joinpath(@__DIR__,"..","..","containers","container-guard.jl"))
    julia = joinpath(Sys.BINDIR,Base.julia_exename())
    marker = "ENTRY_MUST_NOT_EXECUTE"
    for configured in (false,true)
        environment = Dict("PATH"=>ENV["PATH"],"JULIA_DEPOT_PATH"=>join(DEPOT_PATH,':'),
            "JULIA_LOAD_PATH"=>"@:@stdlib","JULIA_NUM_THREADS"=>"1")
        configured && merge!(environment,Dict("LCM_CONTAINER_CPUS"=>"1",
            "LCM_CONTAINER_MEMORY_BYTES"=>string(1024^3),"LCM_CONTAINER_PIDS"=>"128",
            "LCM_CONTAINER_SCRATCH_BYTES"=>string(256*1024^2)))
        output = Pipe(); error_output = Pipe()
        command = setenv(`$julia --startup-file=no --compiled-modules=existing --project=$project --load=$guard -e $("println(\"$marker\")")`,environment)
        child = run(pipeline(command;stdout=output,stderr=error_output);wait=false)
        close(output.in);close(error_output.in)
        out = @async read(output,String)
        err = @async read(error_output,String)
        try
            @test timedwait(()->process_exited(child),30;pollint=0.025) == :ok
            process_exited(child) || kill(child,Base.SIGKILL)
            wait(child)
            @test child.exitcode == 78
            @test !occursin(marker,fetch(out))
            diagnostic = fetch(err)
            @test occursin("LCM container entry denied:",diagnostic)
            @test !occursin("Stacktrace",diagnostic)
            @test !occursin(guard,diagnostic)
            @test ncodeunits(diagnostic) < 256
        finally
            process_exited(child) || kill(child,Base.SIGKILL)
            wait(child)
            close(output);close(error_output);close(child)
            wait(out);wait(err)
        end
    end
end
