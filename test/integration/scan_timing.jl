@testitem "Scan timing / stage ownership, progress, composition and persistence" tags=[:integration, :measurements, :slow] setup=[TestFixtures, SavedResults] begin
    using Logging, Measurements
    const E=LineCableModels.Engine
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0, 1000.0])
    formulation=Formulation()
    plain=compute(problem, formulation)
    fields=(:gridpoint, :wall_seconds, :bytes, :gc_seconds, :compile_seconds, :recompile_seconds)
    timings(log)=[record for record in log.logs if record.group === :timing]
    for timing in (0, 1), progress in (0, 1)

        logger=Test.TestLogger()
        calls=Int[]
        result=with_logger(logger) do
            compute(problem,
                formulation;
                options = (; verbosity = (default = 0, progress, timing),
                    on_result = (_, i,
                        value)->begin
                        @test !haskey(details(value).data, :timing)
                        push!(calls, i)
                    end))
        end
        @test calls == [1]
        @test Z(result) == Z(plain)
        @test Y(result) == Y(plain)
        @test frequencies(result) == frequencies(plain)
        @test !haskey(details(result).data, :timing)
        # A timed point logs one record with its gridpoint and the five measurements.
        @test length(timings(logger)) == timing
        if timing == 1
            record=only(timings(logger)).kwargs
            @test keys(record) == fields
            @test all(>=(0), (record[key] for key in fields[2:end]))
            @test record[:bytes] isa Integer
        end
        logs=filter(record->record.group===:progress, logger.logs)
        @test all(record -> record.group in (:progress, :timing), logger.logs)
        @test isempty(logs) == (progress == 0)
        if progress==1
            @test first(logs).message == "Line parameters computation started"
            @test last(logs).message == "Line parameters computation completed successfully"
            @test last(logs).kwargs[:completed] == 1
            @test last(logs).kwargs[:elapsed_seconds] >= 0
        end
    end
    # Completion attachment retains storage and accommodates other retained fields.
    reused=E.retain_gridpoint(plain, details(plain).data.gridpoint; fields = (reused = true,))
    @test typeof(reused) === typeof(plain)
    @test Z(reused) === Z(plain)
    @test Y(reused) === Y(plain)
    @test details(reused).data.gridpoint === details(plain).data.gridpoint
    for result in (plain, reused)
        @test SavedResults.same(SavedResults.round_trip(result), result)
    end
    timed=(; verbosity = (default = 0, timing = 1))
    calls=Int[]
    forms=Grid((formulation, formulation, formulation))
    batch_timing=Test.TestLogger()
    batch=with_logger(batch_timing) do
        compute(problem, forms; options = (; timed..., on_result = (_, i, _)->push!(calls, i)))
    end
    @test calls == [1, 2, 3]
    @test [record.kwargs[:gridpoint].formulation_index for record in timings(batch_timing)] ==
          [1, 2, 3]
    @test all(value -> Z(value) == Z(plain) && Y(value) == Y(plain), batch)
    @test all(value -> typeof(value) === typeof(first(batch)), batch)
    # Shared-input lookup must not overwrite the throughput clock. A real
    # delayed completion makes an intermediate record available for bounds.
    batch_log=Test.TestLogger()
    with_logger(batch_log) do
        compute(problem, forms; options=(verbosity=(default=0, progress=1),
            on_result=(_, i, _)->(i == 1 && sleep(5.1))))
    end
    intermediate=filter(record->record.message == "Line parameters progress", batch_log.logs)
    @test !isempty(intermediate)
    @test all(record->0 <= record.kwargs[:eta_hours] * 3600 <=
        record.kwargs[:elapsed_seconds] *
        (record.kwargs[:total]-record.kwargs[:completed]), intermediate)
    # Materialization and formulation ordering stay with the existing traversal.
    builds=Ref(0)
    space=Gridspace{LineParametersProblem}(
        temperature->begin
            builds[]+=1
            LineParametersProblem(problem.system; frequencies = problem.frequencies,
                earth_props = problem.earth_props, temperature)
        end,
        (Grid((20.0, 25.0)),))
    choices=Gridspace{typeof(formulation)}(identity, (Grid(Tuple(forms)),))
    log=Test.TestLogger()
    product=with_logger(log) do
        compute(
            ParametricProblem(space,
                ComputationOptions(verbosity = (default = 0, progress = 1, timing = 1))),
            Combinatorial(choices))
    end
    @test builds[] == 2
    @test length(product) == 6
    @test isempty(details(product).data)
    @test length(timings(log)) == 6
    @test all(value -> !haskey(details(value).data, :timing), product)
    # The run writes its lines once, naming the action.
    @test count(record -> occursin("Line parameters computation", string(record.message)),
        log.logs) == 2
    @test last(log.logs).message == "Line parameters computation completed successfully"
    @test last(log.logs).kwargs[:completed] == 6
    uncertain_space=Gridspace{LineParametersProblem}(
        radius->TestFixtures.three_bare_wires_problem(; frequencies = problem.frequencies, radius),
        (Grid((0.04, 0.0425), AbsoluteError(0.0001)),))
    study=ParametricProblem(uncertain_space, ComputationOptions(; timed...))
    linear_log=Test.TestLogger()
    linear=with_logger(() -> compute(study, LinearError(formulation)), linear_log)
    @test isempty(details(linear).data)
    @test length(timings(linear_log)) == 2
    @test SavedResults.same(SavedResults.round_trip(linear), linear)
    mc=MonteCarlo(formulation; trials = 4, seed = 71)
    sampled_log=Test.TestLogger()
    sampled=with_logger(() -> compute(study, mc), sampled_log)
    replay=compute(ParametricProblem(uncertain_space), mc)
    @test statistics(sampled) == statistics(replay)
    @test sampled.point_seeds == replay.point_seeds
    @test sampled.trial_counts == replay.trial_counts == [4, 4]
    # Each trial is a run of one point and logs its own record.
    @test length(timings(sampled_log)) == 8
    @test isempty(details(sampled).data)
    @test samples(sampled) === nothing
    @test histograms(sampled) === nothing
    @test !haskey(details(first(sampled)).data, :timing)
    for retain_details in (false, true), return_samples in (false, true),
        return_histograms in (false, true)
        value=compute(ParametricProblem(uncertain_space),
            MonteCarlo(formulation; trials = 2, seed = 71, retain_details,
                return_samples, return_histograms, bins = 2))
        @test haskey(details(value).data, :trials) == retain_details
        @test SavedResults.same(SavedResults.round_trip(value), value)
    end
    for timing in (0, 1), progress in (0, 1)

        options=ComputationOptions(; verbosity = (default = 0, progress, timing))
        local combination_log=Test.TestLogger()
        values=with_logger(combination_log) do
            (compute(problem, forms; options),
                compute(ParametricProblem(space, options), Combinatorial(choices)),
                compute(ParametricProblem(uncertain_space, options), LinearError(formulation)),
                compute(ParametricProblem(uncertain_space, options), mc))
        end
        @test all(value -> Z(value) == Z(plain) && Y(value) == Y(plain), values[1])
        @test Z.(values[2]) == Z.(product)
        @test Y.(values[2]) == Y.(product)
        @test nominal.(Z(first(values[3]))) == nominal.(Z(first(linear)))
        @test uncertainty.(Z(first(values[3]))) == uncertainty.(Z(first(linear)))
        @test statistics(values[4]) == statistics(sampled)
        @test values[4].point_seeds == sampled.point_seeds
        @test isempty(details(values[4]).data)
        @test length(timings(combination_log)) == timing * (3 + 6 + 2 + 8)
        @test count(record -> endswith(string(record.message), "started"), combination_log.logs) ==
              4progress
        @test count(record -> endswith(string(record.message), "completed successfully"), combination_log.logs) ==
              4progress
    end
    # Rejected physical builds do not log successful timing records.
    rejected=Ref(0)
    rejected_space=Gridspace{LineParametersProblem}(
        radius->begin
            radius>0||(rejected[]+=1; throw(DomainError(radius, "wire radius must be positive")))
            TestFixtures.three_bare_wires_problem(; frequencies = [50.0], radius)
        end,
        (Grid(0.005, AbsoluteError(0.01)),))
    retried_log=Test.TestLogger()
    retried=with_logger(retried_log) do
        compute(ParametricProblem(rejected_space, ComputationOptions(; timed...)),
            MonteCarlo(formulation; trials = 4, seed = 19, retain_details = true,
                max_failures = 50))
    end
    @test rejected[] > 0
    @test length(timings(retried_log)) == only(retried.trial_counts) == 4
    function replace_details(value, retained)
        MonteCarloResult(value.sampling, value.values, value.stats, value.sample_values,
            value.histogram_values, value.root_seed, value.point_seeds, value.trial_counts, retained)
    end
    # A Monte Carlo result retains trials only. Timing is a log of the run.
    @test_throws ArgumentError replace_details(sampled, ComputationDetails(timing = [NamedTuple[], NamedTuple[]]))
    @test_throws ArgumentError replace_details(sampled, ComputationDetails(timing = Vector{NamedTuple}[]))
    @test_throws ArgumentError replace_details(sampled, ComputationDetails(unrelated = [1]))
    for owner in (LineCableModelsCoaxial, LineCableModelsFEM), value in (1, :yes, nothing)

        @test_throws ArgumentError computation_options(owner, ComputationOptions(timing = value))
    end
    for owner in (LineCableModelsCoaxial, LineCableModelsFEM)
        @test !haskey(computation_options(owner, ComputationOptions()).data, :timing)
        @test_throws ArgumentError computation_options(
            owner, ComputationOptions(timing = true, unrelated = true))
    end
    @test_throws ArgumentError computation_options(CableConstantsFormulation, ComputationOptions(timing = true))
    @test_throws ArgumentError compute(ParametricProblem(space, ComputationOptions(timing = 1)), mc)
end
