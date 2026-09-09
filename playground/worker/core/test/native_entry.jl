using Test, UUIDs

@testset "native guard rejects ordinary host entry before scientific code" begin
    project = normpath(joinpath(@__DIR__,".."))
    guard = joinpath(project,"src","native-guard.jl")
    julia = joinpath(Sys.BINDIR,Base.julia_exename())
    marker = "NATIVE_PROFILE_MUST_NOT_LOAD"
    for configured in (false,true)
        environment = Dict("PATH"=>ENV["PATH"],"JULIA_DEPOT_PATH"=>join(DEPOT_PATH,':'),
            "JULIA_LOAD_PATH"=>"@:@stdlib","JULIA_NUM_THREADS"=>"1")
        configured && merge!(environment,Dict("LCM_NATIVE_CPUS"=>"1","LCM_NATIVE_MEMORY_BYTES"=>string(1024^3),
            "LCM_NATIVE_PIDS"=>"128","LCM_NATIVE_SCRATCH_BYTES"=>string(256*1024^2),
            "LCM_NATIVE_UID"=>"1000","LCM_NATIVE_GID"=>"1000",
            "LCM_NATIVE_CGROUP"=>"/user.slice/user-1000.slice/user@1000.service/app.slice/lcm-exec-$(uuid4()).service"))
        output,error_output = Pipe(),Pipe()
        command = setenv(`$julia --startup-file=no --compiled-modules=existing --project=$project --load=$guard -e $("println(\"$marker\")")`,environment)
        child = run(pipeline(command;stdout=output,stderr=error_output);wait=false)
        close(output.in);close(error_output.in)
        out,err = @async(read(output,String)),@async(read(error_output,String))
        try
            @test timedwait(()->process_exited(child),30;pollint=0.025) == :ok
            process_exited(child) || kill(child,Base.SIGKILL)
            wait(child)
            @test child.exitcode == 78
            @test !occursin(marker,fetch(out))
            diagnostic = fetch(err)
            @test occursin("LCM native entry denied:",diagnostic)
            @test !occursin("Stacktrace",diagnostic) && !occursin(guard,diagnostic)
            @test ncodeunits(diagnostic) < 256
        finally
            process_exited(child) || kill(child,Base.SIGKILL)
            wait(child);close(output);close(error_output);close(child);wait(out);wait(err)
        end
    end
end
