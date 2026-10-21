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
    # problem of a scalar problem.
    @test @inferred(compute(problem, Scale{3}())) == CountResult(6)
    @test C.traverse(problem, Scale{3}()) == CountResult(6)
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

    # A parametric problem runs its space with its options, and a Combinatorial
    # formulation sets the retention of details.
    run = compute(C.ParametricProblem(problems),
        C.Combinatorial(Scale{2}(); options = (retain_details = true,)))
    @test collect(run) == CountResult.([2, 4])
    @test length(details(run).data.points) == 2

    # A traversal family member takes its options as keywords, and its positional
    # constructor validates what it builds.
    struct Every{F, O <: C.ComputationOptions} <: C.AbstractTraversal
        inner::F
        options::O
        Every(inner::F, options::O) where {F <: AbstractFormulation, O <: C.ComputationOptions} =
            validate(new{F, O}(inner, options))
    end
    @test Every(Scale{1}(); retain_details = true).options ==
          C.ComputationOptions(; retain_details = true)
    @test_throws ArgumentError Every(Scale{1}(); unused = true)
    @test_throws ArgumentError Every(Scale{1}(), C.ComputationOptions())

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
