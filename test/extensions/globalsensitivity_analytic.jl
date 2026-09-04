@testitem "GlobalSensitivity / Sobol / analytical sensitivity indices" tags=[:extension] begin
    using GlobalSensitivity
    using LineCableModels

    const Grammar = LineCableModels.Grammar
    const PB = LineCableModels.ParametricBuilder

    struct AnalyticProblem <: Grammar.AbstractProblemDefinition
        x1::Float64
        x2::Float64
    end
    struct AnalyticFormulation{Mode} <: Grammar.AbstractFormulation
        solves::Base.RefValue{Int}
    end
    struct AnalyticResult{RValue, LValue}
        resistance::RValue
        inductance::LValue
    end

    LineCableModels.validate(problem::AnalyticProblem) = problem

    function Grammar.compute(
            problem::AnalyticProblem,
            formulation::AnalyticFormulation{Mode};
            options = (;)
    ) where {Mode}
        isempty(options) || error("unexpected analytic options")
        formulation.solves[] += 1
        additive = problem.x1 + problem.x2
        interaction = problem.x1 * problem.x2
        if Mode === :additive
            return AnalyticResult(additive, additive)
        elseif Mode === :interaction
            return AnalyticResult(interaction, interaction)
        elseif Mode === :multioutput
            return AnalyticResult(
                additive,
                reshape([additive, interaction], 2, 1)
            )
        elseif Mode === :normal
            linear = problem.x1 + 2problem.x2
            return AnalyticResult(linear, linear)
        end
        error("unsupported analytic mode")
    end
    function Grammar.observe(result::AnalyticResult, ::typeof(R), indices...)
        return isempty(indices) ? result.resistance : result.resistance[indices...]
    end
    function Grammar.observe(result::AnalyticResult, ::typeof(L), indices...)
        return isempty(indices) ? result.inductance : result.inductance[indices...]
    end
    Grammar.observables(::Type{<:AnalyticResult}) = (R, L)

    function analytic_problem(mu1, sigma1, mu2, sigma2)
        space = PB.Gridspace{AnalyticProblem}(
            AnalyticProblem,
            (
                PB.Grid(mu1, PB.AbsoluteError(sigma1)),
                PB.Grid(mu2, PB.AbsoluteError(sigma2))
            )
        )
        return ParametricProblem(space)
    end

    additive_solves = Ref(0)
    additive = compute(
        analytic_problem(0.0, 1.0, 0.0, 1.0),
        Sensitivity(
            AnalyticFormulation{:additive}(additive_solves),
            GlobalSensitivity.Sobol(),
            (R,);
            samples = 4096,
            distribution = :uniform,
            max_evaluations = 16_384,
            max_output_values = 16_384
        )
    )
    additive_product = only(additive)
    @test additive_solves[] == additive_product.evaluations == 16_384
    @test additive_product.inputs == ["x1", "x2"]
    @test additive_product.active_indices == [1, 2]
    @test only(additive_product.first_order)≈[0.5, 0.5] atol=2.0e-3
    @test only(additive_product.total_order)≈[0.5, 0.5] atol=2.0e-3
    @test additive_product.second_order === nothing
    @test additive_product.confidence == (
        first_order = nothing,
        total_order = nothing,
        second_order = nothing
    )

    interaction_solves = Ref(0)
    interaction = compute(
        analytic_problem(0.0, 1.0, 0.0, 1.0),
        Sensitivity(
            AnalyticFormulation{:interaction}(interaction_solves),
            GlobalSensitivity.Sobol(
                order = [0, 1, 2],
                nboot = 4,
                conf_level = 0.95
            ),
            (R,);
            samples = 4096,
            distribution = :uniform,
            max_evaluations = 98_304,
            max_output_values = 98_304
        )
    )
    interaction_product = only(interaction)
    @test interaction_solves[] == interaction_product.evaluations == 98_304
    @test only(interaction_product.first_order)≈[0.0, 0.0] atol=5.0e-3
    @test only(interaction_product.total_order)≈[1.0, 1.0] atol=5.0e-3
    interaction_second = only(interaction_product.second_order)
    @test interaction_second[1, 2]≈1.0 atol=5.0e-3
    @test interaction_second[2, 1] == 0.0
    @test all(isfinite, only(interaction_product.confidence.first_order))
    @test all(isfinite, only(interaction_product.confidence.total_order))
    @test all(isfinite, only(interaction_product.confidence.second_order))
    @test all(value -> value >= 0, only(interaction_product.confidence.first_order))
    @test all(value -> value >= 0, only(interaction_product.confidence.total_order))
    @test all(value -> value >= 0, only(interaction_product.confidence.second_order))

    multioutput_solves = Ref(0)
    multioutput = compute(
        analytic_problem(0.0, 1.0, 0.0, 1.0),
        Sensitivity(
            AnalyticFormulation{:multioutput}(multioutput_solves),
            GlobalSensitivity.Sobol(order = [0, 1, 2]),
            (@observe(L[:, 1]), R);
            samples = 4096,
            distribution = :uniform,
            max_evaluations = 24_576,
            max_output_values = 73_728
        )
    )
    multioutput_product = only(multioutput)
    @test multioutput_solves[] == multioutput_product.evaluations == 24_576
    @test length(multioutput_product.first_order) == 2
    @test size(multioutput_product.first_order[1]) == (2, 2)
    @test size(multioutput_product.first_order[2]) == (2,)
    @test size(multioutput_product.total_order[1]) == (2, 2)
    @test size(multioutput_product.second_order[1]) == (2, 2, 2)
    @test multioutput_product.first_order[1][1, :]≈[0.5, 0.5] atol=3.0e-3
    @test multioutput_product.first_order[1][2, :]≈[0.0, 0.0] atol=5.0e-3
    @test multioutput_product.total_order[1][2, :]≈[1.0, 1.0] atol=5.0e-3
    @test multioutput_product.first_order[2]≈[0.5, 0.5] atol=3.0e-3

    normal_solves = Ref(0)
    normal = compute(
        analytic_problem(3.0, 1.0, -2.0, 2.0),
        Sensitivity(
            AnalyticFormulation{:normal}(normal_solves),
            GlobalSensitivity.Sobol(),
            (R,);
            samples = 4096,
            distribution = :normal,
            max_evaluations = 16_384,
            max_output_values = 16_384
        )
    )
    normal_product = only(normal)
    expected = [1 / 17, 16 / 17]
    @test normal_solves[] == normal_product.evaluations == 16_384
    @test all(isfinite, only(normal_product.first_order))
    @test all(isfinite, only(normal_product.total_order))
    @test only(normal_product.first_order)≈expected atol=5.0e-3
    @test only(normal_product.total_order)≈expected atol=5.0e-3
end

@testitem "GlobalSensitivity / Sobol / failure boundaries" tags=[:extension] begin
    using GlobalSensitivity
    using LineCableModels
    using QuasiMonteCarlo

    const Grammar = LineCableModels.Grammar
    const PB = LineCableModels.ParametricBuilder

    struct FailureProblem <: Grammar.AbstractProblemDefinition
        x1::Float64
        x2::Float64
    end
    struct FailureFormulation{Mode} <: Grammar.AbstractFormulation
        solves::Base.RefValue{Int}
    end
    struct FailureResult{T}
        value::T
    end
    struct AlternateFailureResult{T}
        value::T
    end

    LineCableModels.validate(problem::FailureProblem) = problem

    function Grammar.compute(
            problem::FailureProblem,
            formulation::FailureFormulation{Mode};
            options = (;)
    ) where {Mode}
        isempty(options) || error("unexpected failure-test options")
        formulation.solves[] += 1
        count = formulation.solves[]
        value = if Mode === :valid
            problem.x1 + problem.x2
        elseif Mode === :complex
            complex(problem.x1, problem.x2)
        elseif Mode === :empty
            Float64[]
        elseif Mode === :nonfinite
            [NaN]
        elseif Mode === :constant
            1.0
        elseif Mode === :shape
            count == 1 ? [problem.x1] : [problem.x1, problem.x2]
        elseif Mode === :type
            problem.x1 + problem.x2
        else
            error("unsupported failure mode")
        end
        Mode === :type && count > 1 && return AlternateFailureResult(value)
        return FailureResult(value)
    end
    Grammar.observe(result::FailureResult, ::typeof(R)) = result.value
    Grammar.observe(result::FailureResult, ::typeof(Z)) = result.value
    Grammar.observe(result::AlternateFailureResult, ::typeof(R)) = result.value
    Grammar.observables(::Type{<:FailureResult}) = (R, Z)
    Grammar.observables(::Type{<:AlternateFailureResult}) = (R,)

    function failure_problem(;
            builder = FailureProblem,
            nominal1 = 0.0,
            sigma1 = 1.0,
            sigma2 = 1.0
    )
        space = PB.Gridspace{FailureProblem}(
            builder,
            (
                PB.Grid(nominal1, PB.AbsoluteError(sigma1)),
                PB.Grid(0.0, PB.AbsoluteError(sigma2))
            )
        )
        return ParametricProblem(space)
    end
    function failure_sensitivity(mode, solves, requests = (R,); kwargs...)
        return Sensitivity(
            FailureFormulation{mode}(solves),
            GlobalSensitivity.Sobol(),
            requests;
            samples = 8,
            distribution = :uniform,
            max_evaluations = 32,
            max_output_values = 32,
            kwargs...
        )
    end
    function capture_exception(f)
        try
            f()
            return nothing
        catch exception
            return exception
        end
    end

    budget_solves = Ref(0)
    budget = failure_sensitivity(
        :valid,
        budget_solves;
        max_evaluations = 63
    )
    error = capture_exception() do
        compute(failure_problem(nominal1 = (0.0, 1.0)), budget)
    end
    @test error isa ArgumentError
    @test occursin("outer points=2", sprint(showerror, error))
    @test occursin("point 1: d=2, evaluations=32", sprint(showerror, error))
    @test occursin("point 2: d=2, evaluations=32", sprint(showerror, error))
    @test occursin("requested total=64", sprint(showerror, error))
    @test budget_solves[] == 0

    scalar_solves = Ref(0)
    scalar = failure_sensitivity(:valid, scalar_solves)
    error = capture_exception() do
        compute(ParametricProblem(FailureProblem(0.0, 0.0)), scalar)
    end
    @test error isa ArgumentError
    @test occursin("problem Gridspace", sprint(showerror, error))
    @test scalar_solves[] == 0

    storage_solves = Ref(0)
    storage = failure_sensitivity(
        :valid,
        storage_solves;
        max_output_values = 31
    )
    error = capture_exception(() -> compute(failure_problem(), storage))
    @test error isa ArgumentError
    @test occursin("after one shape probe", sprint(showerror, error))
    @test storage_solves[] == 1

    inactive_solves = Ref(0)
    inactive = failure_sensitivity(:valid, inactive_solves)
    error = capture_exception() do
        compute(
            failure_problem(sigma1 = 0.0, sigma2 = 0.0),
            inactive
        )
    end
    @test error isa ArgumentError
    @test occursin("use Combinatorial", sprint(showerror, error))
    @test inactive_solves[] == 0

    fixed_solves = Ref(0)
    fixed = failure_sensitivity(
        :valid,
        fixed_solves;
        input_labels = ("active", "fixed")
    )
    fixed_product = only(compute(
        failure_problem(sigma2 = 0.0),
        fixed
    ))
    @test fixed_solves[] == fixed_product.evaluations == 24
    @test fixed_product.inputs == ["active"]
    @test fixed_product.active_indices == [1]
    @test all(isfinite, only(fixed_product.first_order))
    @test all(isfinite, only(fixed_product.total_order))

    labels_solves = Ref(0)
    labels = failure_sensitivity(
        :valid,
        labels_solves;
        input_labels = ("x1",)
    )
    @test_throws DimensionMismatch compute(failure_problem(), labels)
    @test labels_solves[] == 0

    estimator_solves = Ref(0)
    estimator = failure_sensitivity(
        :valid,
        estimator_solves;
        estimator = :Unsupported
    )
    @test_throws ArgumentError compute(failure_problem(), estimator)
    @test estimator_solves[] == 0
    @test_throws ArgumentError failure_sensitivity(
        :valid,
        Ref(0);
        distribution = :lognormal
    )

    for method in (
        GlobalSensitivity.Sobol([0], 1, 0.95),
        GlobalSensitivity.Sobol([0, 1, 3], 1, 0.95),
        GlobalSensitivity.Sobol([0, 1], 0, 0.95),
        GlobalSensitivity.Sobol([0, 1], 1, 0.0),
        GlobalSensitivity.Sobol([0, 1], 1, 1.0),
        GlobalSensitivity.Sobol([0, 1], 1, NaN)
    )
        solves = Ref(0)
        formulation = Sensitivity(
            FailureFormulation{:valid}(solves),
            method,
            (R,);
            samples = 8,
            distribution = :uniform
        )
        @test_throws ArgumentError compute(failure_problem(), formulation)
        @test solves[] == 0
    end

    angle_solves = Ref(0)
    angle = failure_sensitivity(:valid, angle_solves, ((Z, Base.angle),))
    @test_throws ArgumentError compute(failure_problem(), angle)
    @test angle_solves[] == 0

    for (mode, selector, expected_solves, expected_error) in (
        (:complex, Z, 1, ArgumentError),
        (:empty, R, 1, ArgumentError),
        (:nonfinite, R, 1, ArgumentError),
        (:constant, R, 32, ArgumentError),
        (:shape, R, 2, DimensionMismatch),
        (:type, R, 2, ArgumentError)
    )
        solves = Ref(0)
        formulation = failure_sensitivity(mode, solves, (selector,))
        exception = capture_exception(() -> compute(failure_problem(), formulation))
        @test exception isa expected_error
        @test solves[] == expected_solves
    end

    struct RejectPhysicalValues end
    function (::RejectPhysicalValues)(x1, x2)
        abs(x1) <= 0.5 && abs(x2) <= 0.5 || throw(DomainError((x1, x2)))
        return FailureProblem(x1, x2)
    end
    physical_solves = Ref(0)
    physical = failure_sensitivity(:valid, physical_solves)
    @test_throws DomainError compute(
        failure_problem(builder = RejectPhysicalValues()),
        physical
    )
    @test physical_solves[] == 0

    @test_throws ArgumentError PB.UncertainValue(0.0, -1.0)
    @test_throws ArgumentError PB.UncertainValue(Inf, 1.0)
    @test_throws ArgumentError PB.UncertainValue(0.0, Inf)

    unsupported = Sensitivity(
        FailureFormulation{:valid}(Ref(0)),
        GlobalSensitivity.Morris(),
        (R,);
        samples = 8,
        sampler = QuasiMonteCarlo.SobolSample(),
        distribution = :uniform
    )
    @test_throws MethodError compute(failure_problem(), unsupported)
end

@testitem "GlobalSensitivity / Sobol / inference and architecture" tags=[:extension] begin
    using GlobalSensitivity
    using LineCableModels
    using QuasiMonteCarlo

    const Grammar = LineCableModels.Grammar
    extension = Base.get_extension(
        LineCableModels,
        :LineCableModelsGlobalSensitivityExt
    )
    @test extension !== nothing

    struct InferredSensitivityInner <: Grammar.AbstractFormulation end
    formulation = @inferred Sensitivity(
        InferredSensitivityInner(),
        GlobalSensitivity.Sobol(),
        (R,);
        samples = 16
    )
    @test formulation.sampler isa QuasiMonteCarlo.SobolSample
    @test isempty(Test.detect_ambiguities(
        LineCableModels.UQ,
        extension;
        recursive = true
    ))

    @test !isdefined(LineCableModels, :Sobol)
    @test !isdefined(LineCableModels.UQ, :Sobol)
    extension_root = joinpath(
        pkgdir(LineCableModels), "ext", "LineCableModelsGlobalSensitivityExt")
    extension_source = join(
        (read(joinpath(directory, file), String)
        for (directory, _, files) in walkdir(extension_root)
        for file in files if endswith(file, ".jl")),
        '\n'
    )
    @test !occursin(r"struct\s+Sobol\b", extension_source)
    @test !occursin(r"function\s+gsa\b", extension_source)
    @test !occursin("distributed = Val(true)", extension_source)
    @test !occursin(r"ParametricBuilder\._", extension_source)
    @test isempty(filter(
        path -> basename(path) in (
            "helpers.jl",
            "utils.jl",
            "common.jl",
            "manager.jl",
            "provider.jl",
            "context.jl"
        ),
        readdir(extension_root; join = true)
    ))
end
