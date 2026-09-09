@testset "profile loading belongs to explicit work, not bootstrap or inspection" begin
    project = normpath(joinpath(@__DIR__, ".."))
    child = joinpath(@__DIR__, "profile_loading_child.jl")
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    for mode in ("normal", "fail")
        supervisor = ExecutorSupervisor(project; startup_timeout_seconds=30,
            command=`$julia --startup-file=no --compiled-modules=existing --project=$project $child $mode`)
        stages = String[]
        ctx = ExecutionContext("deferred-profile", CancellationToken(),
            (_, stage, _) -> push!(stages, stage), _ -> nothing,
            String[], PreparedResourceCache(), nothing)
        try
            start_executor!(supervisor, ctx)
            @test supervisor.generation == 1
            @test !inspect_preparation!(supervisor, ctx, repeat("0", 64))["ready"]
            @test isempty(stages)
            if mode == "normal"
                prepared = prepare_supervised!(supervisor, ctx, Dict{String,Any}(); timeout_seconds=30)
                @test prepared["evidence"] == Dict("loads"=>1, "value"=>7)
                @test first(stages) == "loading_environment"
                @test inspect_preparation!(supervisor, ctx, prepared["preparation_input_hash"])["ready"]
                repeated = prepare_supervised!(supervisor, ctx, Dict{String,Any}(); timeout_seconds=30)
                @test repeated["cache_status"] == "hit" && repeated["evidence"] == prepared["evidence"]
                spec = OperationSpec("fixture.loaded", identity, (_, _) -> nothing; execution_mode=:supervised)
                @test execute_supervised!(supervisor, spec, ctx, Dict{String,Any}()) ==
                    Dict("loads"=>3, "value"=>7)
            else
                @test_throws RetryableOperationError prepare_supervised!(supervisor, ctx,
                    Dict{String,Any}(); timeout_seconds=30)
                @test supervisor.process === nothing
                @test first(stages) == "loading_environment"
            end
        finally
            stop_executor!(supervisor)
        end
        @test supervisor.process === nothing && isempty(supervisor.io_tasks)
    end
end
