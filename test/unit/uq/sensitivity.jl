@testitem "UQ / Sensitivity / core formulation and result contracts" tags=[:unit] begin
    using LineCableModels

    const Grammar=LineCableModels.Grammar
    const UQ=LineCableModels.UQ

    struct SensitivityInner <: Grammar.AbstractFormulation end
    struct SensitivityMethod end
    struct SensitivitySampler end

    request=@observe(R[1, 1, :])
    formulation=Sensitivity(
        SensitivityInner(),
        SensitivityMethod(),
        (request,);
        samples = 32,
        sampler = SensitivitySampler(),
        input_labels = (:rho, "eps_r"),
        options = (retain_details = true,)
    )
    @test formulation.requests == (request,)
    @test formulation.samples == 32
    @test formulation.distribution === :normal
    @test formulation.estimator === :Jansen1999
    @test formulation.input_labels == ("rho", "eps_r")
    @test formulation.max_evaluations == 200_000
    @test formulation.max_output_values == 20_000_000
    @test formulation.options == (retain_details = true,)
    @test @inferred(computation_options(Sensitivity, (;))) ==
          (retain_details = false,)
    constants=CableConstants(1.0, 2.0, 3.0)
    @test observe(constants, R) == [1.0]

    transformed=Sensitivity(
        SensitivityInner(),
        SensitivityMethod(),
        (@observe((Z, abs)[1, 2, :]),);
        samples = 2,
        sampler = SensitivitySampler(),
        distribution = :uniform,
        estimator = :Sobol2007,
        max_evaluations = 1,
        max_output_values = 1
    )
    @test transformed.distribution === :uniform
    @test transformed.estimator === :Sobol2007
    @test occursin("Sensitivity(samples=2; requests=1)", sprint(show, transformed))
    @test occursin("Global sensitivity attribution", sprint(summary, transformed))

    constructor(args...; kwargs...) = Sensitivity(
        SensitivityInner(),
        SensitivityMethod(),
        args...;
        sampler = SensitivitySampler(),
        kwargs...
    )
    @test_throws ArgumentError constructor((); samples = 2)
    @test_throws ArgumentError constructor(((1, 1),); samples = 2)
    @test_throws ArgumentError constructor(((R, Dict()),); samples = 2)
    @test_throws ArgumentError constructor((R,); samples = 1)
    @test_throws ArgumentError constructor((R,); samples = true)
    @test_throws ArgumentError constructor((R,); samples = 2, distribution = :lognormal)
    @test_throws ArgumentError constructor((R,); samples = 2, estimator = "Jansen1999")
    @test_throws ArgumentError constructor((R,); samples = 2, input_labels = ["x"])
    @test_throws ArgumentError constructor((R,); samples = 2, input_labels = (1,))
    @test_throws ArgumentError constructor((R,); samples = 2, max_evaluations = 0)
    @test_throws ArgumentError constructor((R,); samples = 2, max_output_values = 0)
    @test_throws ArgumentError constructor((R,); samples = 2, options = (unknown = true,))
    @test_throws ArgumentError constructor(
        (R,);
        samples = 2,
        options = (retain_details = :yes,)
    )

    product=(
        inputs = ["rho", "eps_r"],
        active_indices = [1, 2],
        first_order = (reshape([0.6, 0.4], 1, 2),),
        total_order = (reshape([0.7, 0.5], 1, 2),),
        second_order = nothing,
        confidence = (
            first_order = (reshape([0.02, 0.03], 1, 2),),
            total_order = (reshape([0.04, 0.05], 1, 2),),
            second_order = nothing
        ),
        evaluations = 128
    )
    result=SensitivityResult(formulation, [product])
    @test result isa Grammar.AbstractUncertaintyResult{typeof(product)}
    @test Base.IteratorSize(typeof(result)) == Base.HasShape{1}()
    @test Base.IteratorEltype(typeof(result)) == Base.HasEltype()
    @test eltype(typeof(result)) === typeof(product)
    @test length(result) == 1
    @test size(result) == (1,)
    @test firstindex(result) == 1
    @test lastindex(result) == 1
    @test result[1] === product
    @test only(result) === product
    @test collect(result) == [product]
    @test first_order(result) == [product.first_order]
    @test total_order(result) == [product.total_order]
    @test second_order(result) == [nothing]
    @test confidence(result) == [product.confidence]
    @test details(result) == (;)
    @test sprint(show, result) == "SensitivityResult(1 points)"
    @test occursin("Sensitivity result · 1 points", sprint(show, MIME"text/plain"(), result))

    records=[(iterations = 3,)]
    retained=SensitivityResult(
        formulation,
        [product],
        (points = [records],)
    )
    @test details(retained).points[1] === records

    @test_throws ArgumentError SensitivityResult(formulation, NamedTuple[])
    @test_throws MethodError SensitivityResult(formulation, (product,))
    @test_throws DimensionMismatch SensitivityResult(
        formulation,
        [merge(product, (inputs = ["rho"],))]
    )
    @test_throws ArgumentError SensitivityResult(
        formulation,
        [merge(product, (evaluations = 0,))]
    )
    @test_throws DimensionMismatch SensitivityResult(
        formulation,
        [product],
        (points = Vector{NamedTuple}[],)
    )

    @test LineCableModels.Sensitivity === UQ.Sensitivity
    @test LineCableModels.SensitivityResult === UQ.SensitivityResult
    @test LineCableModels.first_order === UQ.first_order
    @test LineCableModels.total_order === UQ.total_order
    @test LineCableModels.second_order === UQ.second_order
end

@testitem "UQ / Sensitivity / extension-unloaded dispatch" tags=[:core_only] begin
    using LineCableModels

    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder

    struct UnloadedInner <: Grammar.AbstractFormulation end
    struct UnloadedMethod end
    struct UnloadedSampler end
    struct UnloadedProblem <: Grammar.AbstractProblemDefinition
        value::Float64
    end

    LineCableModels.validate(problem::UnloadedProblem) = problem

    formulation=Sensitivity(
        UnloadedInner(),
        UnloadedMethod(),
        (@observe(R[1]),);
        samples = 2,
        sampler = UnloadedSampler()
    )
    space=PB.Gridspace{UnloadedProblem}(
        UnloadedProblem,
        (PB.Grid((1.0,)),)
    )
    problem=ParametricProblem(space)
    @test_throws MethodError compute(problem, formulation)
    @test Base.get_extension(
        LineCableModels,
        :LineCableModelsGlobalSensitivityExt
    ) === nothing

    source_root=joinpath(pkgdir(LineCableModels), "src")
    core_source=join(
        (read(joinpath(directory, file), String)
        for (directory, _, files) in walkdir(source_root)
        for file in files if endswith(file, ".jl")),
        '\n'
    )
    @test !occursin("using GlobalSensitivity", core_source)
    @test !occursin("import GlobalSensitivity", core_source)
    @test !occursin("using QuasiMonteCarlo", core_source)
    @test !occursin("import QuasiMonteCarlo", core_source)
end
