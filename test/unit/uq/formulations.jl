@testitem "UQ / study and run options / validated constructors and concrete tuples" tags=[:unit, :uq] begin
    using Random
    inner = CableConstantsFormulation()
    normalize = LineCableModels.formulation_options
    run = LineCableModels.computation_options
    @test LineCableModels.ComputationOptions !== NamedTuple

    # A study is not a formulation. One method normalizes the computation options of the
    # run of every study, and a study's own option is not one of them.
    for owner in (Combinatorial, LinearError, MonteCarlo)
        @test !(owner(inner) isa AbstractFormulation)
        @test run(owner, ComputationOptions()) == ComputationOptions(;
            verbosity = (default = 0,), on_result = nothing, retain_details = false)
        @test run(owner, ComputationOptions(retain_details = true)).data.retain_details
        @test run(owner, ComputationOptions(verbosity = (default = 1,))).data.verbosity ==
              (default = 1,)
        @test_throws ArgumentError run(owner, ComputationOptions(unused = true))
        @test_throws ArgumentError run(owner, ComputationOptions(trials = 8))
        @test_throws ArgumentError run(owner, ComputationOptions(retain_details = 1))
    end
    # Combinatorial and LinearError have no options of their own.
    for owner in (Combinatorial, LinearError)
        @test @inferred(owner(inner)).inner === inner
        @test_throws MethodError owner(inner; retain_details = true)
        @test_throws MethodError owner(inner, ComputationOptions())
    end
    @test_throws MethodError Combinatorial(Grid((1, 2)))
    @test_throws MethodError Combinatorial(Gridspace{Int}(identity, (Grid((1, 2)),)))

    # The keyword constructor of MonteCarlo merges its options, and `formulation_options` of
    # its type normalizes them. A computation option of the run is not one of them.
    keyword = @inferred MonteCarlo(inner; options=(trials=8,))
    @test @inferred(MonteCarlo(inner; trials=8)) === keyword
    @test @inferred(MonteCarlo(inner, keyword.options)) === keyword
    @test normalize(MonteCarlo, keyword.options) === keyword.options
    @test normalize(keyword) === keyword.options
    @test_throws ArgumentError MonteCarlo(inner; options=(unused=true,))
    @test_throws ArgumentError MonteCarlo(inner, FormulationOptions(unused=true))
    @test_throws ArgumentError MonteCarlo(inner; retain_details=true)
    @test_throws MethodError MonteCarlo(inner, ComputationOptions(trials=8))
    # A `Grid` among the options makes a space of option records, which MonteCarlo refuses.
    @test_throws MethodError MonteCarlo(inner; trials=Grid((4, 8)))
    @test_throws MethodError MonteCarlo(inner; options=(seed=Grid((1, 2)),))
    @test_throws MethodError MonteCarlo(inner, FormulationOptions(distribution=Grid((:normal, :uniform))))
    # Explicit type parameters must not restore the generated raw constructor.
    invalid = FormulationOptions(unused=true)
    @test_throws MethodError MonteCarlo{typeof(inner), typeof(invalid)}(inner, invalid)
    @test_throws TypeError MonteCarlo{typeof(inner), NamedTuple}

    plain = (trials=Int16(8), confidence=0.9f0, cdf_tol=0.1f0,
        distribution=:uniform, seed=Int16(42), return_samples=true,
        return_histograms=true, bins=Int16(2), max_failures=Int16(5))
    normalized = @inferred normalize(MonteCarlo, FormulationOptions(plain))
    @test normalized.data.trials === 8
    @test normalized.data.seed === UInt64(42)
    @test normalized.data.confidence === Float64(0.9f0)
    @test normalized.data.cdf_tol === Float64(0.1f0)
    @test normalized.data.bins === 2
    @test normalized.data.max_failures === 5
    @test (@inferred MonteCarlo(inner; options=plain)).options === normalized
    @test (@inferred MonteCarlo(inner; plain...)).options === normalized
    @test (@inferred MonteCarlo(inner; trials=8, options=(seed=42,))).options ===
          (@inferred MonteCarlo(inner; options=(trials=8, seed=42))).options
    @test_throws ArgumentError MonteCarlo(inner; trials=8, options=(trials=8,))
    @test_throws ArgumentError MonteCarlo(inner; unused=true)
    @test_throws ArgumentError normalize(MonteCarlo, FormulationOptions((distribution=:unsupported,)))

    for key in (:trials, :bins, :max_failures)
        for value in (true, 0, -1, 1.5, "2", big(typemax(Int)) + 1)
            options = NamedTuple{(key,)}((value,))
            @test_throws ArgumentError normalize(MonteCarlo, FormulationOptions(options))
        end
    end
    @test_throws ArgumentError normalize(MonteCarlo, FormulationOptions((max_failures=nothing,)))
    for key in (:confidence, :cdf_tol)
        for value in (true, 0, 1, -0.5, NaN, Inf, "0.5",
                big(1) - eps(BigFloat), big(2)^(-2000))
            @test_throws ArgumentError normalize(MonteCarlo, FormulationOptions(NamedTuple{(key,)}((value,))))
        end
    end
    for key in (:return_samples, :return_histograms)
        @test_throws ArgumentError normalize(MonteCarlo, FormulationOptions(NamedTuple{(key,)}((1,))))
    end
    for seed in (true, -1, 1.5, "42", big(typemax(UInt64)) + 1)
        @test_throws ArgumentError normalize(MonteCarlo, FormulationOptions((; seed)))
    end
    @test normalize(MonteCarlo, FormulationOptions((seed=0,))).data.seed === UInt64(0)
    @test normalize(MonteCarlo, FormulationOptions((seed=typemax(UInt64),))).data.seed === typemax(UInt64)
    automatic = @inferred MonteCarlo(inner)
    @test automatic.options.data.trials === nothing
    @test automatic.options.data.seed === nothing
    @test automatic.options.data.bins === nothing

    sampler = (rng, center, deviation) -> center + deviation * randn(rng)
    custom = @inferred MonteCarlo(inner; options=(distribution=sampler,))
    @test custom.options.data.distribution === sampler
    @test fieldtype(typeof(custom.options.data), :distribution) === typeof(sampler)
end

@testitem "UQ / studies / a study is not a point, and options follow two levels" tags=[:unit, :uq] setup=[TestFixtures] begin
    using Measurements, Logging
    frequencies = [50.0]
    problem = TestFixtures.three_bare_wires_problem(; frequencies)
    space = Gridspace{LineParametersProblem}(
        rho -> TestFixtures.three_bare_wires_problem(; rho, frequencies),
        (Grid((10.0, 100.0)),))
    f = Formulation()
    mc = MonteCarlo(f; trials = 4, seed = 1)

    # Grids vary the inputs of one scalar problem and formulation. The calls that process
    # spaces take none, a study contains a core formulation and a parametric problem
    # contains core problems: each other form is refused at its call.
    @test_throws MethodError MonteCarlo(f; trials = Grid((4, 8)))
    @test_throws MethodError Combinatorial(f; retain_details = Grid((true, false)))
    @test_throws MethodError Combinatorial(mc)
    @test_throws MethodError Combinatorial(mc; retain_details = true)
    @test_throws MethodError MonteCarlo(Combinatorial(f))
    @test_throws MethodError LinearError(mc)
    @test_throws MethodError MonteCarlo(LinearError(f))
    @test_throws MethodError Combinatorial(Grid((mc, MonteCarlo(f; trials = 8))))
    @test_throws MethodError ParametricProblem(ParametricProblem(space))
    @test_throws MethodError ParametricProblem(
        Grid((ParametricProblem(space), ParametricProblem(space))))
    @test_throws MethodError compute(space, Grid((mc, MonteCarlo(f; trials = 8))))
    @test_throws MethodError compute(ParametricProblem(space), Combinatorial(mc))
    @test_throws MethodError compute(ParametricProblem(space), LinearError(mc))
    @test_throws MethodError compute(ParametricProblem(space),
        MonteCarlo(LinearError(f); trials = 4, seed = 1))
    @test_throws MethodError compute(ParametricProblem(ParametricProblem(space)), f)
    @test_throws MethodError compute(
        Grid((ParametricProblem(space), ParametricProblem(space))), f)
    @test_throws ArgumentError compute(space, f;
        options = (verbosity = Grid(((default = 0,), (default = 1,))),))
    formulations = Grid((f, Formulation()))
    @test_throws MethodError MonteCarlo(formulations)
    @test_throws MethodError MonteCarlo(formulations; trials = 4)
    @test_throws MethodError LinearError(formulations)
    @test_throws MethodError MonteCarlo(formulations, ComputationOptions())
    @test_throws MethodError Combinatorial(f,
        Grid((ComputationOptions(), ComputationOptions())))
    @test_throws MethodError ParametricProblem(Grid((space, space)))

    # The options of the backend for every point are those of the parametric problem, and
    # the computation options of the run those of the call: each reaches its owner.
    traced = ParametricProblem(space, ComputationOptions(trace = true))
    log = Test.TestLogger(min_level = Logging.Info)
    result = with_logger(log) do
        compute(traced, Combinatorial(f);
            options = (verbosity = (default = 0, progress = 1), retain_details = true))
    end
    @test all(value -> haskey(details(value).data, :trace), result)
    @test length(details(result).data.points) == 2
    @test count(record -> record.group === :progress &&
        endswith(string(record.message), "computation started"), log.logs) == 1
    @test count(record -> record.group === :progress &&
        endswith(string(record.message), "completed successfully"), log.logs) == 1
    untraced = compute(ParametricProblem(space), Combinatorial(f))
    @test !any(value -> haskey(details(value).data, :trace), untraced)
    @test isempty(details(untraced).data)
    # Each kind of option is refused at the other level, and no backend takes the options of
    # the run, which a scalar call gives to its run too.
    @test_throws ArgumentError compute(traced, Combinatorial(f); options = (trace = true,))
    @test_throws ArgumentError compute(ParametricProblem(space,
        ComputationOptions(retain_details = true)), Combinatorial(f))
    @test_throws ArgumentError compute(ParametricProblem(space,
        ComputationOptions(verbosity = (default = 1,))), Combinatorial(f))
    calls = Int[]
    constants = compute(CableConstantsProblem(TestFixtures.coaxial_design()),
        CableConstantsFormulation();
        options = (verbosity = (default = 0,), on_result = (_, index, _) -> push!(calls, index)))
    @test constants isa CableConstants
    @test calls == [1]

    # `compute(space, x; options)` gives each key to its owner, for a core formulation, whose
    # default study is `Combinatorial`, and for a study alike.
    split = compute(space, f; options = (trace = true, retain_details = true))
    @test split isa ParametricResult
    @test all(value -> haskey(details(value).data, :trace), split)
    @test length(details(split).data.points) == 2
    uncertain = Gridspace{LineParametersProblem}(
        rho -> TestFixtures.three_bare_wires_problem(; rho, frequencies),
        (Grid(100.0, AbsoluteError(5.0)),))
    calls = Int[]
    sampled = compute(uncertain, mc;
        options = (retain_details = true, on_result = (_, index, _) -> push!(calls, index)))
    @test sampled isa MonteCarloResult
    @test calls == [1]
    @test length.(details(sampled).data.trials) == [4]
    @test sampled.point_seeds ==
          compute(ParametricProblem(uncertain), mc).point_seeds == [UInt64(1)]
    linear = compute(uncertain, LinearError(f); options = (retain_details = true,))
    @test linear isa LinearErrorResult
    @test length(details(linear).data.points) == 1
end

@testitem "UQ / studies / a study chained into modal analysis is the chain inside each point" tags=[:unit, :uq] setup=[TestFixtures] begin
    using Measurements
    frequencies = [50.0, 5000.0]
    f = Formulation()
    modal = ModalAnalysisFormulation(:default)
    same(a, b) = isequal(Z(a), Z(b)) && isequal(Y(a), Y(b)) && isequal(gamma(a), gamma(b))
    pointwise(problems) = [compute(ModalAnalysisProblem(compute(p, f)), modal) for p in problems]

    # A Combinatorial or a LinearError study gives each point's line parameters to modal
    # analysis: the same values as the chain of each point.
    space = Gridspace{LineParametersProblem}(
        rho -> TestFixtures.three_bare_wires_problem(; rho, frequencies),
        (Grid((10.0, 100.0)),))
    expected = pointwise(collect(space))
    for chained in (compute(space, Combinatorial(f); modal = :default),
            compute(space, Combinatorial(f), modal),
            compute(ParametricProblem(space), Combinatorial(f), modal))
        @test all(same.(collect(chained), expected))
    end
    uncertain = Gridspace{LineParametersProblem}(
        rho -> TestFixtures.three_bare_wires_problem(; rho, frequencies),
        (Grid((10.0, 100.0), AbsoluteError(5.0)),))
    expected = pointwise([LineCableModels.materialize(point)
                          for point in LineCableModels.points(uncertain)])
    for chained in (compute(uncertain, LinearError(f); modal = :default),
            compute(uncertain, LinearError(f), modal),
            compute(ParametricProblem(uncertain), LinearError(f), modal))
        @test all(same.(collect(chained), expected))
    end

    # A Monte Carlo study chained into modal analysis is refused at the call, because Monte
    # Carlo of the modal outputs is not built.
    mc = MonteCarlo(f; trials = 2, seed = 1)
    for refused in (() -> compute(uncertain, mc; modal = :default),
            () -> compute(uncertain, mc, modal),
            () -> compute(ParametricProblem(uncertain), mc, modal))
        caught = try
            refused()
            nothing
        catch exception
            exception
        end
        @test caught isa ArgumentError
        @test occursin("Monte Carlo of the modal outputs is not built", sprint(showerror, caught))
    end
end
