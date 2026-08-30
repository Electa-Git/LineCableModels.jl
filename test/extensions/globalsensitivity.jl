@testitem "Extensions / GlobalSensitivity / dependency API contract" tags=[:extension] begin
    using GlobalSensitivity
    using QuasiMonteCarlo

    @test pkgversion(GlobalSensitivity) == v"2.12.1"
    @test pkgversion(QuasiMonteCarlo) == v"0.3.6"

    samples = 16
    dimension = 2
    bootstraps = 2
    matrices = QuasiMonteCarlo.generate_design_matrices(
        samples,
        dimension,
        SobolSample(),
        2bootstraps
    )
    @test length(matrices) == 2bootstraps
    @test all(matrix -> size(matrix) == (dimension, samples), matrices)
    @test all(matrix -> all(isfinite, matrix), matrices)
    @test all(matrix -> all(value -> 0 <= value <= 1, matrix), matrices)

    A = reduce(hcat, matrices[1:bootstraps])
    B = reduce(hcat, matrices[(bootstraps + 1):(2bootstraps)])
    callback_shapes = Tuple{Int, Int}[]
    function batch_model(points)
        push!(callback_shapes, size(points))
        return vcat(
            (points[1, :] .+ points[2, :])',
            (points[1, :] .* points[2, :])'
        )
    end

    first_total = GlobalSensitivity.gsa(
        batch_model,
        GlobalSensitivity.Sobol(order = [0, 1], nboot = bootstraps),
        A,
        B;
        batch = true,
        Ei_estimator = :Jansen1999,
        distributed = Val(false)
    )
    @test callback_shapes == [(dimension, samples * bootstraps * (dimension + 2))]
    @test size(first_total.S1) == (2, dimension)
    @test size(first_total.ST) == (2, dimension)
    @test size(first_total.S1_Conf_Int) == (2, dimension)
    @test size(first_total.ST_Conf_Int) == (2, dimension)
    @test first_total.S2 === nothing
    @test first_total.S2_Conf_Int === nothing

    empty!(callback_shapes)
    second = GlobalSensitivity.gsa(
        batch_model,
        GlobalSensitivity.Sobol(order = [0, 1, 2], nboot = bootstraps),
        A,
        B;
        batch = true,
        Ei_estimator = :Jansen1999,
        distributed = Val(false)
    )
    @test callback_shapes == [(dimension, samples * bootstraps * (2dimension + 2))]
    @test size(second.S1) == (2, dimension)
    @test size(second.ST) == (2, dimension)
    @test size(second.S2) == (dimension, dimension, 2)
    @test size(second.S1_Conf_Int) == (2, dimension)
    @test size(second.ST_Conf_Int) == (2, dimension)
    @test size(second.S2_Conf_Int) == (dimension, dimension, 2)

    for estimator in (:Homma1996, :Sobol2007, :Jansen1999, :Janon2014)
        estimate = GlobalSensitivity.gsa(
            batch_model,
            GlobalSensitivity.Sobol(),
            matrices[1],
            matrices[2];
            batch = true,
            Ei_estimator = estimator,
            distributed = Val(false)
        )
        @test size(estimate.S1) == (2, dimension)
        @test size(estimate.ST) == (2, dimension)
    end
end

@testitem "Extensions / GlobalSensitivity / Sobol integration contract" tags=[:extension] begin
    using Distributions
    using GlobalSensitivity
    using LineCableModels
    using QuasiMonteCarlo

    const Grammar = LineCableModels.Grammar
    const PB = LineCableModels.ParametricBuilder

    struct SobolContractProblem <: Grammar.AbstractProblemDefinition
        x1::Float64
        x2::Float64
    end
    struct SobolContractFormulation <: Grammar.AbstractFormulation end
    struct SobolContractResult
        outputs::Vector{Float64}
    end

    const contract_solves = Ref(0)
    function Grammar.compute(
            problem::SobolContractProblem,
            ::SobolContractFormulation;
            options = (;)
    )
        isempty(options) || error("unexpected core options")
        contract_solves[] += 1
        return SobolContractResult([
            problem.x1 + problem.x2,
            problem.x1 * problem.x2
        ])
    end
    Grammar.observe(result::SobolContractResult, ::typeof(R)) = result.outputs
    function Grammar.computation_details(
            ::Val{SobolContractFormulation},
            result::SobolContractResult
    )
        return (checksum = sum(result.outputs),)
    end

    space = PB.Gridspace{SobolContractProblem}(
        SobolContractProblem,
        (
            PB.Grid(0.0, PB.AbsoluteError(1.0)),
            PB.Grid(0.0, PB.AbsoluteError(1.0))
        )
    )
    problem = ParametricProblem(space)
    method = GlobalSensitivity.Sobol(order = [0, 1], nboot = 2)
    formulation = Sensitivity(
        SobolContractFormulation(),
        method,
        (R,);
        samples = 32,
        distribution = :uniform,
        input_labels = (:left, :right),
        max_evaluations = 256,
        max_output_values = 512,
        options = (retain_details = true,)
    )
    @test formulation.sampler isa QuasiMonteCarlo.SobolSample
    @test Base.get_extension(
        LineCableModels,
        :LineCableModelsGlobalSensitivityExt
    ) !== nothing

    result = compute(problem, formulation)
    @test result isa SensitivityResult
    @test contract_solves[] == 256
    @test length(result) == 1
    product = only(result)
    @test product.inputs == ["left", "right"]
    @test product.active_indices == [1, 2]
    @test product.evaluations == 256
    @test length(product.first_order) == 1
    @test size(only(product.first_order)) == (2, 2)
    @test size(only(product.total_order)) == (2, 2)
    @test product.second_order === nothing
    @test size(only(product.confidence.first_order)) == (2, 2)
    @test size(only(product.confidence.total_order)) == (2, 2)
    @test product.confidence.second_order === nothing
    @test length(details(result).evaluations) == 1
    @test length(only(details(result).evaluations)) == product.evaluations
    @test eltype(only(details(result).evaluations)) <: NamedTuple
end
