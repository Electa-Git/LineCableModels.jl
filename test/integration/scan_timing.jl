@testitem "Scan timing / stage ownership, progress, composition and persistence" tags=[:integration] setup=[TestFixtures] begin
    using Logging, Measurements, JSON3, Serialization
    const E=LineCableModels.Engine
    const IE=LineCableModels.ImportExport
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0, 1000.0])
    formulation=Formulation()
    plain=compute(problem, formulation)
    fields=(:wall_seconds, :bytes, :gc_seconds, :compile_seconds, :recompile_seconds)
    timed_result=Ref{Any}(nothing)
    for timing in (false, true), progress in (0, 1)

        logger=Test.TestLogger()
        calls=Int[]
        result=with_logger(logger) do
            compute(problem,
                formulation;
                options = (; timing, verbosity = (default = 0, progress),
                    on_result = (_, i,
                        value)->begin
                        @test haskey(details(value).data, :timing) == timing
                        push!(calls, i)
                    end))
        end
        @test calls == [1]
        @test Z(result) == Z(plain)
        @test Y(result) == Y(plain)
        @test frequencies(result) == frequencies(plain)
        @test haskey(details(result).data, :timing) == timing
        if timing
            timed_result[]=result
            record=details(result).data.timing
            @test keys(record) == fields
            @test all(>=(0), values(record))
            @test record.bytes isa Integer
        end
        logs=filter(record->record.group===:progress, logger.logs)
        @test all(record -> record.group === :progress, logger.logs)
        @test isempty(logs) == (progress == 0)
        if progress==1
            @test first(logs).message == "Line parameters computation started"
            @test last(logs).message == "Line parameters computation completed successfully"
            @test last(logs).kwargs[:completed] == 1
            @test last(logs).kwargs[:elapsed_seconds] >= 0
        end
    end
    timed=timed_result[]
    # Completion attachment retains storage and accommodates empty reuse records.
    reused=E.retain_gridpoint(timed, details(timed).data.gridpoint; fields = (timing = (;),))
    @test typeof(reused) === typeof(timed)
    @test Z(reused) === Z(timed)
    @test Y(reused) === Y(timed)
    @test details(reused).data.gridpoint === details(timed).data.gridpoint
    for result in (plain, timed, reused)
        restored=IE.deserialize_value(JSON3.read(JSON3.write(IE.serialize_value(result)), Dict{
            String, Any}))
        @test Z(restored) == Z(result)
        @test get(details(restored).data, :timing, nothing) ==
              get(details(result).data, :timing, nothing)
    end
    calls=Int[]
    forms=[formulation, formulation, formulation]
    batch=compute(problem, forms; options = (
        timing = true, on_result = (_, i, _)->push!(calls, i)))
    @test calls == [1, 2, 3]
    @test all(value -> keys(details(value).data.timing) == fields, batch)
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
                ComputationOptions(timing = true, verbosity = (default = 0, progress = 1))),
            Combinatorial(choices))
    end
    @test builds[] == 2
    @test length(product) == 6
    @test isempty(details(product).data)
    @test all(value -> haskey(details(value).data, :timing), product)
    @test !any(record -> occursin("Line parameters computation", string(record.message)), log.logs)
    @test last(log.logs).message == "Parametric computation completed successfully"
    @test last(log.logs).kwargs[:completed] == 6
    uncertain_space=Gridspace{LineParametersProblem}(
        radius->TestFixtures.three_bare_wires_problem(; frequencies = problem.frequencies, radius),
        (Grid((0.04, 0.0425), AbsoluteError(0.0001)),))
    study=ParametricProblem(uncertain_space, ComputationOptions(timing = true))
    linear=compute(study, LinearError(formulation))
    @test isempty(details(linear).data)
    @test haskey(details(first(linear)).data, :timing)
    restored=IE.deserialize_value(JSON3.read(JSON3.write(IE.serialize_value(linear)), Dict{
        String, Any}))
    @test details(first(restored)).data.timing == details(first(linear)).data.timing
    @test nominal.(Z(first(restored))) == nominal.(Z(first(linear)))
    @test uncertainty.(Z(first(restored))) == uncertainty.(Z(first(linear)))
    mc=MonteCarlo(formulation; trials = 4, seed = 71)
    sampled=compute(study, mc)
    replay=compute(ParametricProblem(uncertain_space), mc)
    @test statistics(sampled) == statistics(replay)
    @test sampled.point_seeds == replay.point_seeds
    @test sampled.trial_counts == replay.trial_counts == [4, 4]
    @test keys(details(sampled).data) == (:timing,)
    @test length.(details(sampled).data.timing) == [4, 4]
    @test samples(sampled) === nothing
    @test histograms(sampled) === nothing
    @test !haskey(details(first(sampled)).data, :timing)
    for retain_details in (false, true), return_samples in (false, true),
        return_histograms in (false, true)
        value=compute(study,
            MonteCarlo(formulation; trials = 2, seed = 71, retain_details,
                return_samples, return_histograms, bins = 2))
        local restored_mc=IE.deserialize_value(JSON3.read(
            JSON3.write(IE.serialize_value(value)), Dict{String, Any}))
        @test details(restored_mc).data.timing == details(value).data.timing
        @test statistics(restored_mc) == statistics(value)
        @test haskey(details(restored_mc).data, :trials) == retain_details
    end
    # The scientific record codec is also valid through Julia serialization.
    buffer=IOBuffer()
    serialize(buffer, IE.serialize_value(sampled))
    seekstart(buffer)
    @test details(IE.deserialize_value(deserialize(buffer))).data.timing ==
          details(sampled).data.timing
    for timing in (false, true), progress in (0, 1)

        options=ComputationOptions(; timing, verbosity = (default = 0, progress))
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
        @test haskey(details(values[4]).data, :timing) == timing
        @test count(record -> endswith(string(record.message), "started"), combination_log.logs) ==
              4progress
        @test count(record -> endswith(string(record.message), "completed successfully"), combination_log.logs) ==
              4progress
    end
    # Rejected physical builds do not create successful timing entries.
    rejected_space=Gridspace{LineParametersProblem}(
        radius->begin
            radius>0||throw(DomainError(radius, "wire radius must be positive"))
            TestFixtures.three_bare_wires_problem(; frequencies = [50.0], radius)
        end,
        (Grid(0.005, AbsoluteError(0.01)),))
    retried=compute(ParametricProblem(rejected_space, ComputationOptions(timing = true)),
        MonteCarlo(formulation; trials = 4, seed = 19, retain_details = true,
            on_error = :retry, max_failures = 50))
    @test !isempty(only(details(retried).data.failures))
    @test length(only(details(retried).data.timing)) == only(retried.trial_counts) == 4
    function replace_details(value, retained)
        MonteCarloResult(value.formulation, value.values, value.stats, value.sample_values,
            value.histogram_values, value.root_seed, value.point_seeds, value.trial_counts, retained)
    end
    @test_throws DimensionMismatch replace_details(sampled, ComputationDetails(timing = [NamedTuple[], NamedTuple[]]))
    @test_throws DimensionMismatch replace_details(sampled, ComputationDetails(timing = Vector{NamedTuple}[]))
    @test_throws ArgumentError replace_details(sampled, ComputationDetails(unrelated = [1]))
    for owner in (LineCableModelsCoaxial, LineCableModelsFEM), value in (1, :yes, nothing)

        @test_throws ArgumentError computation_options(owner, ComputationOptions(timing = value))
    end
    for owner in (LineCableModelsCoaxial, LineCableModelsFEM)
        @test !computation_options(owner, ComputationOptions()).data.timing
        @test_throws ArgumentError computation_options(
            owner, ComputationOptions(timing = true, unrelated = true))
    end
    @test_throws ArgumentError computation_options(CableConstantsFormulation, ComputationOptions(timing = true))
    @test_throws ArgumentError compute(ParametricProblem(space, ComputationOptions(timing = 1)), mc)
end

@testitem "Scan timing / benchmark repeats completed projections" tags=[:integration] setup=[TestFixtures] begin
    using Measurements
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0])
    formulation=Formulation()
    timed=compute(problem, formulation; options = (timing = true,))
    reused=LineCableModels.Engine.retain_gridpoint(
        timed, details(timed).data.gridpoint; fields = (timing = (;),))
    calls=Ref(0)
    record=LineCableModels.benchmark(; samples = 3, warmup = 2) do
        calls[]+=1
        isodd(calls[]) ? timed : reused
    end
    @test calls[] == 5
    @test record.result === timed
    @test record.timings == [details(timed).data.timing, (;), details(timed).data.timing]
    builds = Ref(0)
    space = Gridspace{LineParametersProblem}(radius -> begin
        builds[] += 1
        TestFixtures.three_bare_wires_problem(; frequencies=[50.], radius)
    end, (Grid((0.04, 0.0425)),))
    study = ParametricProblem(space, ComputationOptions(timing=true))
    for result in ([timed, reused],
        compute(study, Combinatorial(Grid((formulation, formulation, formulation)))),
        compute(study, LinearError(formulation)),
        compute(study, MonteCarlo(formulation; trials = 2, seed = 1)))
        before = builds[]
        projected=LineCableModels.benchmark(()->result; samples = 2)
        @test builds[] == before
        @test projected.result === result
        @test projected.timings[1] == projected.timings[2]
        @test projected.timings[1] !== projected.timings[2]
        if result isa MonteCarloResult
            @test projected.timings[1] == details(result).data.timing
            @test projected.timings[1] !== details(result).data.timing
            @test projected.timings[1][1] !== details(result).data.timing[1]
        else
            @test projected.timings[1] == [details(core).data.timing for core in result]
        end
    end
    for count in (0, -1, true, 1.5)
        @test_throws ArgumentError LineCableModels.benchmark(() -> timed; samples = count)
    end
    for count in (-1, true, 1.5)
        @test_throws ArgumentError LineCableModels.benchmark(() -> timed; warmup = count)
    end
    @test_throws r"timing=true" LineCableModels.benchmark(() -> compute(problem, formulation))
    @test_throws r"timing=true" LineCableModels.benchmark(() -> [
        timed, compute(problem, formulation)])
    calls[]=0
    @test_throws ErrorException LineCableModels.benchmark(; samples = 3) do
        calls[] += 1
        error("failure")
    end
    @test calls[] == 1
end
