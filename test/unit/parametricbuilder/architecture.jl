@testitem "ParametricBuilder / public API / explicit materialized boundary" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    @test LineCableModels.Material === PB.Material
    @test LineCableModels.Material === LineCableModels.Materials.Material
    @test LineCableModels.MaterialsLibrary === LineCableModels.Materials.MaterialsLibrary
    @test LineCableModels.CablesLibrary === LineCableModels.DataModel.CablesLibrary
    @test LineCableModels.Fortescue === LineCableModels.Engine.Transforms.Fortescue
    @test LineCableModels.Conductor === PB.Conductor
    @test !isdefined(LineCableModels, :Tubular)
    @test !isdefined(LineCableModels, :ConductorGroup)
    @test !isdefined(LineCableModels, :FormulationSet)
    @test !isdefined(LineCableModels, :Thickness)
    @test !isdefined(LineCableModels, :Diameter)
    @test isdefined(LineCableModels.DataModel, :Tubular)
    @test isdefined(LineCableModels.Materials, :Material)
end

@testitem "Materials / public API / library insertion" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    library=LineCableModels.MaterialsLibrary(add_defaults = false)
    copper=LineCableModels.Material(
        rho = 1.7241e-8,
        eps_r = 1.0,
        mu_r = 0.999994,
        T0 = 20.0,
        alpha = 0.00393
    )

    @test copper isa LineCableModels.Material
    add!(library, :copper, copper)

    @test library["copper"] isa LineCableModels.Materials.Material
    @test library["copper"].rho == 1.7241e-8
    @test_throws ArgumentError LineCableModels.Material(rho = 1.0, combine = :outer)
    @test_throws ArgumentError add!(
        library,
        :sweep,
        LineCableModels.Material(rho = Grid((1.0, 2.0)))
    )
    @test_throws ArgumentError add!(
        library,
        :uncertain,
        LineCableModels.Material(rho = Grid(1.0, 1.0))
    )
end

@testitem "ParametricBuilder / CableBuilder / numeric radius and thickness staging" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    design_spec=PB.CableBuilder(
        "typed-cable",
        PB.Conductor.Solid(:core; radius = Grid((0.010, 0.012)), material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe),
        nominal = (;
            designation_code = "GRID-CABLE",
            U0 = 18.0,
            U = 30.0,
            conductor_cross_section = 1000.0
        )
    )

    designs=collect(design_spec)
    @test length(designs) == 2
    @test all(design -> design isa LineCableModels.DataModel.CableDesign, designs)
    @test designs[1].components[1].conductor_group.r_ex == 0.010
    @test designs[1].nominal_data.designation_code == "GRID-CABLE"
    @test designs[1].nominal_data.conductor_cross_section == 1000.0
    @test designs[1].components[1].insulator_group.r_in == 0.010
    @test designs[1].components[2].conductor_group.r_in == 0.014
    @test designs[1].components[2].insulator_group.r_in == 0.015

    earth=PB.Earth(
        rho = Grid((10.0, 100.0, 1000.0)),
        eps_r = Grid((100.0, 10.0, 5.0));
        combine = :zip
    )
    @test [(item.rho, item.eps_r) for item in earth] == [
        (10.0, 100.0),
        (100.0, 10.0),
        (1000.0, 5.0)
    ]

    frequencies_grid=Grid(([50.0], [50.0, 500.0]))
    systems=PB.SystemBuilder(
        "typed-system",
        first(designs),
        PB.trifoil(y = -1.0, spacing = 0.05, phases = :core=>(1, 2, 3));
        earth = PB.Earth(rho = 100.0),
        frequencies = frequencies_grid
    )
    @test length(systems) == 2
    @test map(problem -> problem.frequencies, collect(systems)) ==
          [[50.0], [50.0, 500.0]]
    first_system=first(systems).system
    first_position, second_position=first_system.cables[1:2]
    @test hypot(
        first_position.horz - second_position.horz,
        first_position.vert - second_position.vert
    ) ≈ 0.05
end

@testitem "Computation / compute! / deterministic parametric and Monte Carlo ownership" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Random
    using Distributions
    import LineCableModels.ParametricBuilder as PB
    import LineCableModels.Computation as Computation

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    uncertain_radius=Grid((0.010, 0.012), 2.0; key = :core_radius)
    design_spec=PB.CableBuilder(
        "compute-cable",
        PB.Conductor.Solid(:core; radius = uncertain_radius, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )
    problem=CableConstantsProblem(design_spec)
    formulation=Formulation()

    @test_throws ArgumentError compute!(problem, formulation)

    direct=compute!(problem, formulation; run = FullParametric())
    @test direct isa FullParametricResult{<:CableConstants}
    @test length(direct) == 2
    @test all(value -> value.R isa Measurement, direct)

    policy=MonteCarlo(
        trials = 12,
        seed = 0x1234,
        return_samples = true,
        return_histograms = true
    )
    monte_carlo=compute!(problem, formulation; run = policy)
    @test monte_carlo isa FullParametricResult{<:MonteCarloResult{<:CableConstants}}
    @test length(monte_carlo) == 2
    @test all(value -> value.trials == 12, monte_carlo)
    @test all(value -> value.samples !== nothing, monte_carlo)
    @test all(value -> value.histograms !== nothing, monte_carlo)
    @test all(value -> value.uncertain isa PB.UncertainValue, monte_carlo)
    @test length(manifest(monte_carlo).hash) == 64
    @test manifest(monte_carlo).execution_policy.actual_root_seed == UInt64(0x1234)

    constants_frame=DataFrame(first(monte_carlo))
    @test constants_frame.quantity == ["R", "L", "C"]
    @test constants_frame.trials == fill(12, 3)
    @test constants_frame.cdf_tol == fill(policy.cdf_tol, 3)
    @test !(:ci_half in propertynames(constants_frame))
    @test DataFrames.metadata(constants_frame, "monte_carlo").manifest_hash ==
          manifest(first(monte_carlo)).hash

    resistance_pdf=first(monte_carlo).histograms.R
    @test cdf(resistance_pdf, maximum(resistance_pdf)) == 1.0
    @test pdf(resistance_pdf, quantile(resistance_pdf, 0.5)) >= 0
    pdf_draw=rand(MersenneTwister(8), resistance_pdf)
    @test first(resistance_pdf.edges) <= pdf_draw <= last(resistance_pdf.edges)

    repeated=compute!(problem, formulation; run = policy)
    @test manifest(repeated).hash == manifest(monte_carlo).hash
    @test map(value -> value.representation, repeated) ==
          map(value -> value.representation, monte_carlo)

    distribution_problem=CableConstantsProblem(PB.CableBuilder(
        "distribution-cable",
        PB.Conductor.Solid(
            :core;
            radius = Grid(0.010, AbsoluteError(1e-4)),
            material = copper
        ),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    ))
    distribution_run=compute!(
        distribution_problem,
        formulation;
        run = MonteCarlo(trials = 4, seed = 7, distribution = Laplace())
    )
    @test distribution_run isa MonteCarloResult{<:CableConstants}
    @test distribution_run.seed == UInt64(7)
    @test manifest(distribution_run).execution_policy.actual_root_seed == UInt64(7)
    descriptor=only(Grid(1.0, AbsoluteError(0.1)))
    @test_throws ArgumentError rand(
        MersenneTwister(9),
        descriptor;
        distribution = Cauchy()
    )
    automatic_run=compute!(
        distribution_problem,
        formulation;
        run = MonteCarlo(
            trials = nothing,
            confidence = 0.5,
            cdf_tol = 0.9,
            seed = 8
        )
    )
    @test automatic_run.trials == Computation._dkw_trials(3, 0.5, 0.9)

    correlated=Measurements.measurement(first(monte_carlo))
    @test correlated isa CableConstants{<:Measurement}
    @test Measurements.cov(correlated.R, correlated.L) != 0

    @test Computation._dkw_trials(3, 0.95, 0.02) ==
          ceil(Int, log(2 * 3 / 0.05) / (2 * 0.02^2))
end

@testitem "Engine / Formulation / materialized and Gridspace line problems" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    formulation=Formulation()

    deterministic_design=PB.CableBuilder(
        "line-cable",
        PB.Conductor.Solid(
            :core;
            radius = Grid((0.010, 0.011)),
            material = copper
        ),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )
    deterministic_problem=PB.SystemBuilder(
        "line-system",
        deterministic_design,
        PB.at(x = 0.0, y = -1.0, phases = (:core=>1, :screen=>0));
        earth = PB.Earth(rho = 100.0),
        frequencies = [50.0]
    )

    ordinary=compute!(first(deterministic_problem), formulation)
    @test ordinary isa LineParameters
    parametric=compute!(
        deterministic_problem,
        formulation;
        options = (verbosity = 0,)
    )
    @test parametric isa FullParametricResult{<:LineParameters}
    @test length(parametric) == 2

    uncertain_design=PB.CableBuilder(
        "uncertain-line-cable",
        PB.Conductor.Solid(
            :core;
            radius = Grid(0.010, 1.0),
            material = copper
        ),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )
    uncertain_problem=PB.SystemBuilder(
        "uncertain-line-system",
        uncertain_design,
        PB.at(x = 0.0, y = -1.0, phases = (:core=>1, :screen=>0));
        earth = PB.Earth(rho = 100.0),
        frequencies = [50.0]
    )

    @test_throws ArgumentError compute!(uncertain_problem, formulation)
    direct=compute!(uncertain_problem, formulation; run = FullParametric())
    @test direct isa FullParametricResult{<:LineParameters}
    @test only(direct).Z.values[1] isa Complex{<:Measurement}

    monte_carlo=compute!(
        uncertain_problem,
        formulation;
        run = MonteCarlo(
            trials = 3,
            seed = 11,
            return_samples = true,
            return_histograms = true,
            bins = 2
        )
    )
    @test monte_carlo isa MonteCarloResult{<:LineParameters}
    @test monte_carlo.trials == 3
    @test statistics(monte_carlo) isa RLCG
    @test samples(monte_carlo) isa RLCG
    @test histograms(monte_carlo) isa RLCG
    @test first(histograms(monte_carlo).R) isa HistogramPDF
    @test Measurements.measurement(monte_carlo) isa LineParameters

    line_frames=DataFrame(monte_carlo)
    @test size(line_frames) == size(statistics(monte_carlo).R)
    @test line_frames[1, 1, 1].quantity == ["R", "L", "C", "G"]
    @test all(==("Ω/km"), line_frames[1, 1, 1].unit[1:1])
end

@testitem "Computation / invalid configurations / strict and auditable handling" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    design=PB.CableBuilder(
        "invalid-config-cable",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )
    problem=CableConstantsProblem(design; separation = Grid((-1.0, -2.0)))

    @test_throws ArgumentError compute!(problem, Formulation(); run = FullParametric())
    skipped=compute!(
        problem,
        Formulation();
        run = FullParametric(invalid = :skip)
    )
    @test skipped isa FullParametricResult
    @test isempty(skipped)
    @test length(skipped.failures) == 2
    @test all(failure -> failure.exception_type == "ArgumentError", skipped.failures)
end

@testitem "DataModel / DataFrame / presentation does not compute" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    design=only(PB.CableBuilder(
        "presentation-cable",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    ))

    @test names(DataFrame(design))[1] == "property"
    components=DataFrame(design, :components)
    detailed=DataFrame(design, :detailed)
    @test components.property == [
        :radius_in_con,
        :radius_ext_con,
        :rho_con,
        :alpha_con,
        :mu_con,
        :radius_ext_ins,
        :eps_ins,
        :mu_ins,
        :loss_factor_ins
    ]
    @test ncol(components) == length(design.components) + 1
    @test detailed.property[1:6] ==
          ["type", "r_in", "r_ex", "diam_in", "diam_ext", "thickness"]
    @test ncol(detailed) ==
          1 + sum(
        length(component.conductor_group.layers) +
        length(component.insulator_group.layers)
    for component in design.components
    )
    @test all(
        column -> all(!ismissing, column[2:6]),
        eachcol(detailed[:, Not(:property)])
    )
    @test_throws ArgumentError DataFrame(design, :baseparams)
    @test_throws ErrorException DataFrame(design, :unsupported)
    constants=compute!(CableConstantsProblem(design), Formulation())
    rendered=DataFrame(constants)
    @test rendered.parameter == ["R", "L", "C"]
    @test rendered.value == [constants.R, constants.L, constants.C]
end
