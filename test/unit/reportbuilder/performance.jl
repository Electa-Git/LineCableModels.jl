@testitem "ReportBuilder / recorded performance has explicit scopes and no collection" tags=[:unit] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition, tabulate
    using DataFrames
    definition=BenchmarkTableDefinition()
    @test all(isempty, values(tabulate(definition, nothing)))
    source=(backend = :pscad, scope = :compile_call, seconds = 0.0328,
        excludes = (:output_readiness, :transfer), reused = false)
    execution=(
        reference = (backend = :pscad,
            timing = (scope = :compute_call_wall, seconds = 12.0, source_timings = source),
            reused = true, execution_wall_seconds = 0.01),
        candidate = (backend = :fem,
            timing = (scope = :compute_call_wall, seconds = 4.0,
                source_timings = (backend = :getdp, scope = :worker_sum,
                    solve_seconds = 7.0, recovered_columns = 2))))
    measure=(scope = :compute_call_wall, median_seconds = 10.0, bytes = 0, samples = 1,
        observations = [(
            seconds = 10.0, bytes = 0, reused = false, source_timings = source)],
        environment = (threads = 1, workers = 2, instrumented = false),
        calculation = (input = "reference workload", trials = 512),
        policy = (progress = false, callbacks = false, allocation_scope = :julia,
            allocation_statistic = :maximum, warmup = :owned_call))
    performance=(reference = measure,
        candidate = merge(
            measure, (median_seconds = 1.0,
                calculation = (input = "candidate workload", trials = 0))),
        speedup = 10.0, comparable = false, passes = nothing,
        settings = (minimum_speedup = 2.0, samples = 3, seconds = 20.0))
    tables=tabulate(definition,
        (execution = execution, performance = performance,
            checksum_verified=missing,workload_verified=true,session=(id="original",)))
    @test tables.execution.reused[1]===true
    @test tables.execution.execution_wall_seconds[1]==0.01
    @test tables.source_timings.scope==[:compile_call, :worker_sum]
    @test tables.source_timings.seconds[1]==0.0328
    @test tables.source_timings.solve_seconds[2]==7.0
    @test tables.performance.samples==[1, 1]
    @test all(iszero, tables.performance.allocated_bytes)
    @test all(==(:maximum), tables.performance.allocation_statistic)
    @test all(ismissing,tables.performance.checksum_verified)
    @test all(==(true),tables.performance.workload_verified)
    @test tables.performance.session_id==["original","original"]
    @test only(tables.performance_comparison.reference_over_candidate)==10.0
    @test !only(tables.performance_comparison.comparable)
    @test ismissing(only(tables.performance_comparison.passes))
    @test only(tables.performance_comparison.requested_samples)==3
    @test size(tables.performance_samples, 1)==2
    # Counts retain their measurement meaning: timed repetitions, not MC trials.
    @test tables.performance.requested_samples==[3,3]
    for table in values(tables),column in eachcol(table),value in column
        @test value isa Union{Number,Bool,Symbol,AbstractString,Missing}
    end
end
