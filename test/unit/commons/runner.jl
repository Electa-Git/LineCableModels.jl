@testitem "Commons / runner / one point, a space and a parametric problem" tags=[:unit, :commons] begin
    using Logging
    const C = LineCableModels.Commons
    struct CountProblem <: AbstractProblemDefinition
        value::Int
    end
    LineCableModels.validate(problem::CountProblem) = problem
    struct Scale{K} <: AbstractFormulation end
    LineCableModels.computation_options(::Type{<:Scale}, options::C.ComputationOptions) = options
    struct CountResult <: AbstractCoreResult
        value::Int
    end
    function LineCableModels.compute(resolved::C.Resolved{<:CountProblem}, ::Scale{K};
            options = C.ComputationOptions()) where {K}
        return CountResult(K * resolved.problem.value)
    end
    LineCableModels.computation_details(::Type{<:Scale}, ::CountResult) = C.ComputationDetails()

    # The default resolve gives each formulation `Resolved(problem, (;))`, a tuple for a tuple.
    problem = CountProblem(2)
    options = C.ComputationOptions()
    @test C.resolve(problem, (Scale{1}(), Scale{3}()); options) ==
          (C.Resolved(problem, (;)), C.Resolved(problem, (;)))
    @test C.resolve(problem, [Scale{1}()]; options) isa Vector{<:C.Resolved}

    # A scalar call is a run of one point and returns its value, as does the parametric
    # problem of a scalar problem with a core formulation.
    @test @inferred(compute(problem, Scale{3}())) == CountResult(6)
    @test @inferred(C.traverse(C.ParametricProblem(problem), Scale{3}())) == CountResult(6)
    @test @inferred(compute(C.ParametricProblem(problem), Scale{3}())) == CountResult(6)
    one = compute(problem, C.Grid((Scale{1}(), Scale{3}())))
    @test collect(one) == CountResult.([2, 6])
    @test one[1, 2] == CountResult(6)

    # A space returns its values problem index fastest. `on_result` runs in computation
    # order, problem by problem, with each point's index in the space.
    indices = Int[]
    problems = C.Gridspace{CountProblem}(CountProblem, (C.Grid((1, 2)),))
    result = compute(problems, C.Grid((Scale{1}(), Scale{10}()));
        options = (on_result = (_, index, _) -> push!(indices, index),))
    @test result isa C.ParametricResult
    @test collect(result) == CountResult.([1, 2, 10, 20])
    @test indices == [1, 3, 2, 4]
    @test result[2, 2] == CountResult(20)
    @test result.gridpoints == CartesianIndices((2, 2))

    # A parametric problem runs its space with its options, and the computation options of
    # the run, in the call, set the retention of details.
    run = compute(C.ParametricProblem(problems), C.Combinatorial(Scale{2}());
        options = (retain_details = true,))
    @test collect(run) == CountResult.([2, 4])
    @test length(details(run).data.points) == 2

    # A study of one's own holds its inner formulation, and the computation options of its
    # run are normalized by the one method of every study. A study is not a formulation,
    # and a parametric problem is a problem space.
    struct Every{F} <: C.AbstractTraversal
        inner::F
        Every(inner::F) where {F <: AbstractFormulation} = validate(new{F}(inner))
    end
    @test computation_options(Every, C.ComputationOptions(retain_details = true)) ==
          C.ComputationOptions(; verbosity = (default = 0,), on_result = nothing,
        retain_details = true)
    @test_throws ArgumentError computation_options(Every, C.ComputationOptions(unused = true))
    @test !(Every(Scale{1}()) isa AbstractFormulation)
    @test C.ParametricProblem(problem) isa C.AbstractProblemSpace
    @test !(C.ParametricProblem(problem) isa AbstractProblemDefinition)
    @test_throws MethodError C.ParametricProblem(C.ParametricProblem(problem))
    @test_throws MethodError C.Combinatorial(C.Combinatorial(Scale{1}()))

    # At level 1 in the `timing` group, each point logs one record with its gridpoint and
    # the five timing fields, and the result does not record timing.
    timing_log = Test.TestLogger(min_level = Logging.Info)
    timed = Logging.with_logger(timing_log) do
        compute(problems, C.Grid((Scale{1}(), Scale{10}()));
            options = (verbosity = (default = 0, timing = 1),))
    end
    @test collect(timed) == CountResult.([1, 2, 10, 20])
    records = [record for record in timing_log.logs if record.group === :timing]
    @test [(record.kwargs[:gridpoint].problem_index, record.kwargs[:gridpoint].formulation_index)
           for record in records] == [(1, 1), (1, 2), (2, 1), (2, 2)]
    @test all(record -> keys(record.kwargs) == (:gridpoint, :wall_seconds, :bytes,
        :gc_seconds, :compile_seconds, :recompile_seconds), records)
    # A formulation type of one's own defines its normalizer, here one that keeps the options.
    @test computation_options(Scale{1}, options) === options

    # The default returns a result type without retained details as it is.
    @test C.retain_gridpoint(CountResult(1), C.gridpoint_id()) == CountResult(1)
end

@testitem "Commons / runner / a space records a failed point and continues" tags=[:unit, :commons] begin
    const C = LineCableModels.Commons
    struct FailProblem <: AbstractProblemDefinition
        value::Int
    end
    LineCableModels.validate(problem::FailProblem) = problem
    struct OtherProblem <: AbstractProblemDefinition
        value::Int
    end
    LineCableModels.validate(problem::OtherProblem) = problem
    struct Weight{K} <: AbstractFormulation end
    LineCableModels.computation_options(::Type{<:Weight}, options::C.ComputationOptions) =
        options
    struct WeightResult <: AbstractCoreResult
        value::Int
    end
    struct Refused <: Exception end
    # The point of problem 3 is outside the domain of the heavier formulation, and problem 5
    # refuses every formulation with an exception of its own.
    function LineCableModels.compute(resolved::C.Resolved{<:FailProblem}, ::Weight{K};
            options = C.ComputationOptions()) where {K}
        value = resolved.problem.value
        value == 3 && K == 10 && throw(DomainError(value, "outside the domain"))
        value == 5 && throw(Refused())
        return WeightResult(K * value)
    end
    LineCableModels.compute(resolved::C.Resolved{<:OtherProblem}, ::Weight{K};
        options = C.ComputationOptions()) where {K} = WeightResult(100K * resolved.problem.value)
    LineCableModels.computation_details(::Type{<:Weight}, result::WeightResult) =
        C.ComputationDetails((value = result.value,))
    # The cables of problem 2 overlap: its build throws.
    overlapping = value -> value == 2 ?
        throw(OverlapError((1, 2), -1e-3, "cable exteriors overlap")) : FailProblem(value)
    space = C.Gridspace{FailProblem}(overlapping, (C.Grid((1, 2, 3, 4)),))
    formulations = C.Grid((Weight{1}(), Weight{10}()))

    # The study records the overlapping problem for each formulation and the point outside
    # the domain, in the order of the space. It computes every other point, and `on_result`
    # runs for the computed points only.
    indices = Int[]
    result = compute(space, formulations;
        options = (on_result = (_, index, _) -> push!(indices, index),))
    @test [(gridpoint.problem_index, gridpoint.formulation_index)
           for (gridpoint, _) in result.failures] == [(2, 1), (2, 2), (3, 2)]
    @test [typeof(exception) for (_, exception) in result.failures] ==
          [OverlapError, OverlapError, DomainError]
    @test result.failures isa Vector{Pair{typeof(C.gridpoint_id()), Exception}}
    @test indices == [1, 5, 3, 4, 8]
    @test length(result) + length(result.failures) == length(result.gridpoints) == 8
    # The values equal those of a run without the failed points. A scalar point infers.
    alone = [@inferred(compute(FailProblem(value), weight))
             for (value, weight) in ((1, Weight{1}()), (3, Weight{1}()), (4, Weight{1}()),
        (1, Weight{10}()), (4, Weight{10}()))]
    @test collect(result) == alone
    @test [result[1, 1], result[3, 1], result[4, 1], result[1, 2], result[4, 2]] == alone
    lookup = try
        result[2, 1]
    catch exception
        exception
    end
    @test lookup isa ErrorException
    @test occursin("problem 2 and formulation 1", sprint(showerror, lookup))
    @test occursin("cable exteriors overlap", sprint(showerror, lookup))
    @test_throws ErrorException result[3, 2]
    @test sprint(show, result) == "ParametricResult(5 values, 3 failed)"
    @test occursin("5 values, 3 failed", sprint(show, MIME("text/plain"), result))
    # Retained details cover the computed points.
    retained = compute(C.ParametricProblem(space), C.Combinatorial(formulations);
        options = (retain_details = true,))
    @test [record.data.value for record in details(retained).data.points] ==
          [1, 3, 4, 10, 40]

    # A space whose points build two problem types computes as a run of each, and records
    # the failed points among them.
    mixed = value -> value == 2 ? OtherProblem(value) : overlapping(value)
    both = compute(C.Gridspace{AbstractProblemDefinition}(mixed, (C.Grid((1, 2, 4)),)),
        formulations)
    @test isempty(both.failures)
    @test collect(both) == WeightResult.([1, 200, 4, 10, 2000, 40])
    recorded = compute(C.Gridspace{AbstractProblemDefinition}(
            value -> value == 4 ? OtherProblem(value) : overlapping(value),
            (C.Grid((1, 2, 3, 4)),)),
        formulations)
    @test [(gridpoint.problem_index, gridpoint.formulation_index)
           for (gridpoint, _) in recorded.failures] == [(2, 1), (2, 2), (3, 2)]
    @test collect(recorded) == WeightResult.([1, 3, 400, 10, 4000])

    # An exception that the study does not record stops it with its own backtrace.
    struct Broken <: AbstractFormulation end
    LineCableModels.computation_options(::Type{<:Broken}, options::C.ComputationOptions) =
        options
    broken_root(value) = sqrt("$(value)")
    LineCableModels.compute(resolved::C.Resolved{<:FailProblem}, ::Broken;
        options = C.ComputationOptions()) =
        resolved.problem.value == 4 ? broken_root(4) : WeightResult(0)
    stopped, trace = try
        compute(C.Gridspace{FailProblem}(FailProblem, (C.Grid((1, 4)),)), Broken())
    catch exception
        exception, catch_backtrace()
    end
    @test stopped isa MethodError
    @test any(frame -> frame.func === :broken_root, stacktrace(trace))

    # A space with no computed point throws, with the count and the first failure.
    every = try
        compute(C.Gridspace{FailProblem}(overlapping, (C.Grid((2, 2)),)), Weight{1}())
    catch exception
        exception
    end
    @test every isa ErrorException
    @test occursin("all 2 points", sprint(showerror, every))
    @test occursin("cable exteriors overlap", sprint(showerror, every))

    # A user records an exception of their own with one method.
    refusing = C.Gridspace{FailProblem}(FailProblem, (C.Grid((1, 5)),))
    @test_throws Refused compute(refusing, Weight{1}())
    LineCableModels.on_error(::C.Combinatorial, exception::Refused) = exception
    own = compute(refusing, Weight{1}())
    @test collect(own) == [WeightResult(1)]
    @test last(only(own.failures)) === Refused()

    # The owners of the exception types declare what a study records. Monte Carlo draws
    # again after an infeasible draw, and any other exception stops a study.
    study = C.Combinatorial(Weight{1}())
    domain = DomainError(1, "outside the domain")
    @test @inferred(C.on_error(study, domain)) === domain
    overlap = OverlapError((1, 2), -1e-3, "overlap")
    @test @inferred(C.on_error(study, overlap)) === overlap
    @test_throws ArgumentError try
        throw(ArgumentError("not recorded"))
    catch exception
        C.on_error(study, exception)
    end
    # Linear propagation records what Combinatorial records, the user's exception included.
    linear = LinearError(Formulation())
    @test C.on_error(linear, overlap) === overlap
    @test C.on_error(linear, Refused()) === Refused()
    sampled = MonteCarlo(Formulation())
    @test C.on_error(sampled, overlap) === overlap
    @test C.on_error(sampled, domain) === domain
    @test_throws Refused try
        throw(Refused())
    catch exception
        C.on_error(sampled, exception)
    end
end
