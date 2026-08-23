@testitem "ParametricBuilder / public API / explicit materialized boundary" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine
    const Computation=LineCableModels.Computation

    @test LineCableModels.AbstractProblemDefinition === Grammar.AbstractProblemDefinition
    @test LineCableModels.AbstractFormulation === Grammar.AbstractFormulation
    @test LineCableModels.AbstractProblemResult === Grammar.AbstractProblemResult
    @test LineCableModels.AbstractParametricResult === Grammar.AbstractParametricResult
    @test LineCableModels.AbstractUncertaintyResult === Grammar.AbstractUncertaintyResult
    @test LineCableModels.compute === Grammar.compute === Engine.compute
    @test LineCableModels.observables === Grammar.observables
    @test LineCableModels.primitives === Grammar.primitives
    @test LineCableModels.preprocess === Grammar.preprocess
    @test LineCableModels.FormulationOptions === Grammar.FormulationOptions === NamedTuple
    @test LineCableModels.ComputationOptions === Grammar.ComputationOptions === NamedTuple
    for name in (
        :Grid, :AbsoluteError, :DeterministicGrid, :RelativeGrid, :AbsoluteGrid,
        :AbstractGrid, :AbstractUncertainGrid, :UncertainValue,
        :Gridspace, :Configuration
    )
        @test getproperty(LineCableModels, name) === getproperty(Grammar, name)
        @test getproperty(PB, name) === getproperty(Grammar, name)
    end
    for name in (
        :Combinatorial, :LinearError, :MonteCarlo, :ParametricProblem,
        :ParametricResult, :LinearErrorResult, :MonteCarloResult,
        :CalculationManifest, :ConfigurationFailure, :SampleSummary,
        :HistogramDensity, :RLCG
    )
        @test getproperty(LineCableModels, name) === getproperty(Grammar, name)
        @test getproperty(Computation, name) === getproperty(Grammar, name)
    end
    @test getproperty(PB, Symbol("@gridspace")) ===
          getproperty(Grammar, Symbol("@gridspace"))
    @test getproperty(PB, Symbol("@relax")) ===
          getproperty(Grammar, Symbol("@relax"))
    @test Computation.Combinatorial === Grammar.Combinatorial
    @test Computation.ParametricResult === Grammar.ParametricResult
    @test Engine.AbstractProblemDefinition === Grammar.AbstractProblemDefinition
    @test Engine.AbstractFormulation === Grammar.AbstractFormulation

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
    for removed in (
        :ProblemDefinition, :EMTFormulation, :EMTTrace, :ComputeOptions,
        :AbstractRunPolicy, :FullParametric, :FullParametricResult,
        :AbstractSpec, :MaterialSpec, :PartSpec, :CableDesignSpec,
        :PositionSpec, :EarthSpec, :SystemSpec, :HistogramPDF
    )
        @test !isdefined(LineCableModels, removed)
    end
    @test PB.MaterialDefinition !== nothing
    @test PB.PartDefinition !== nothing
    @test PB.CableDesignDefinition !== nothing
    @test PB.PositionDefinition !== nothing
    @test PB.EarthDefinition !== nothing
    @test PB.SystemDefinition !== nothing
    @test isdefined(LineCableModels.DataModel, :Tubular)
    @test isdefined(LineCableModels.Materials, :Material)
end

@testitem "Engine / grammar / strict analytical and computation options" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine
    struct TestFormulationFamily<:Engine._FormulationFamily end

    formulation=Formulation()
    @test formulation isa AnalyticalFormulation
    @test formulation.backend == Val(:analytical)
    @test formulation.options.output == Val(:parameters)
    @test Formulation(:analytical; options = (output = :trace,)).options.output ==
          Val(:trace)
    @test_throws MethodError Formulation(:EMT)
    @test_throws MethodError Formulation(:unknown)
    @test_throws ArgumentError Formulation(:analytical; options = (unknown = true,))
    @test_throws ArgumentError Formulation(:analytical; options = (output = :unknown,))
    @test_throws ArgumentError Formulation(:analytical; options = (kron_reduction = :yes,))

    family=Engine._LineParameterFormulationFamily
    @test Engine._active_formulation_backend(family) === :analytical
    @test Engine._activate_formulation_backend!(TestFormulationFamily, :unimplemented) ===
          :unimplemented
    @test Engine._active_formulation_backend(TestFormulationFamily) === :unimplemented
    @test Engine._active_formulation_backend(family) === :analytical
    @test Formulation(:analytical) isa AnalyticalFormulation
    @test Formulation() isa AnalyticalFormulation

    options=Engine.computation_options((
        verbosity = (default = 1, NLsolve = 0),
        output_basis = :total
    ))
    @test options isa ComputationOptions
    @test options.output_basis == Val(:total)
    @test Engine.verbosity(options, :NLsolve) == 0
    @test Engine.verbosity(options, :unlisted) == 1
    @test_throws ArgumentError Engine.computation_options((unknown = true,))
    @test_throws ArgumentError Engine.computation_options((output_basis = :unknown,))
    @test isempty(methods(Grammar.preprocess))
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

    designs=collect(Gridspace(design_spec))
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
    @test [(item.layers[end].rho, item.layers[end].eps_r) for item in Gridspace(earth)] == [
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
    system_space=Gridspace(systems)
    @test length(system_space) == 2
    @test map(problem -> problem.frequencies, collect(system_space)) ==
          [[50.0], [50.0, 500.0]]
    first_system=first(system_space).system
    first_position, second_position=first_system.cables[1:2]
    @test hypot(
        first_position.horz - second_position.horz,
        first_position.vert - second_position.vert
    ) ≈ 0.05
end

@testitem "Computation / compute / deterministic parametric and Monte Carlo ownership" tags=[:unit] setup=[
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
    space=CableConstantsProblem(design_spec)
    formulation=Formulation()

    @test_throws MethodError compute(space, formulation)

    direct=compute(ParametricProblem(space), Combinatorial(formulation))
    @test direct isa ParametricResult{<:CableConstants}
    @test length(direct) == 2
    @test all(value -> value.R isa Measurement, direct)

    fixed_design=PB.CableBuilder(
        "fixed-compute-cable",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )
    singleton=compute(
        ParametricProblem(CableConstantsProblem(fixed_design)),
        Combinatorial(formulation)
    )
    @test singleton isa ParametricResult{<:CableConstants}
    @test length(singleton) == 1

    policy=MonteCarlo(
        formulation;
        trials = 12,
        seed = 0x1234,
        return_samples = true,
        return_histograms = true
    )
    monte_carlo=compute(ParametricProblem(space), policy)
    @test monte_carlo isa MonteCarloResult{<:CableConstants}
    @test length(monte_carlo) == 2
    @test monte_carlo.details[:random].trials == [12, 12]
    @test length(samples(monte_carlo)) == 2
    @test length(histograms(monte_carlo)) == 2
    @test length(manifest(monte_carlo).hash) == 64
    @test monte_carlo.details[:random].root_seed == UInt64(0x1234)

    resistance_pdf=first(histograms(monte_carlo)).R
    @test cdf(resistance_pdf, maximum(resistance_pdf)) == 1.0
    @test pdf(resistance_pdf, quantile(resistance_pdf, 0.5)) >= 0
    pdf_draw=rand(MersenneTwister(8), resistance_pdf)
    @test first(resistance_pdf.edges) <= pdf_draw <= last(resistance_pdf.edges)

    repeated=compute(ParametricProblem(space), policy)
    @test manifest(repeated).hash == manifest(monte_carlo).hash
    @test result(repeated) == result(monte_carlo)

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
    distribution_run=compute(
        ParametricProblem(distribution_problem),
        MonteCarlo(formulation; trials = 4, seed = 7, distribution = Laplace(),
            return_samples = true, return_histograms = true)
    )
    @test distribution_run isa MonteCarloResult{<:CableConstants}
    @test distribution_run.details[:random].root_seed == UInt64(7)
    constants_frame=DataFrame(distribution_run)
    @test constants_frame.quantity == ["R", "L", "C"]
    @test constants_frame.trials == fill(4, 3)
    @test constants_frame.cdf_tol == fill(0.02, 3)
    @test !(:ci_half in propertynames(constants_frame))
    @test DataFrames.metadata(constants_frame, "monte_carlo").manifest_hash ==
          manifest(distribution_run).hash
    descriptor=only(Grid(1.0, AbsoluteError(0.1)))
    @test_throws ArgumentError rand(
        MersenneTwister(9),
        descriptor;
        distribution = Cauchy()
    )
    automatic_run=compute(
        ParametricProblem(distribution_problem),
        MonteCarlo(
            formulation;
            trials = nothing,
            confidence = 0.5,
            cdf_tol = 0.9,
            seed = 8
        )
    )
    @test only(automatic_run.details[:random].trials) ==
          Computation._dkw_trials(3, 0.5, 0.9)

    propagated=compute(ParametricProblem(space), LinearError(formulation))
    @test propagated isa LinearErrorResult{<:CableConstants}
    @test Measurements.cov(first(propagated).R, first(propagated).L) != 0
    @test !applicable(Measurements.measurement, monte_carlo)

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

    deterministic_space=Gridspace(deterministic_problem)
    ordinary=compute(first(deterministic_space), formulation)
    @test ordinary isa LineParameters
    parametric=compute(
        ParametricProblem(
            deterministic_space,
            (verbosity = (default = 0,),)
        ),
        Combinatorial(formulation)
    )
    @test parametric isa ParametricResult{<:LineParameters}
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

    uncertain_space=Gridspace(uncertain_problem)
    @test_throws MethodError compute(uncertain_space, formulation)
    direct=compute(ParametricProblem(uncertain_space), LinearError(formulation))
    @test direct isa LinearErrorResult{<:LineParameters}
    @test only(direct).Z.values[1] isa Complex{<:Measurement}

    monte_carlo=compute(
        ParametricProblem(uncertain_space),
        MonteCarlo(
            formulation;
            trials = 3,
            seed = 11,
            return_samples = true,
            return_histograms = true,
            bins = 2
        )
    )
    @test monte_carlo isa MonteCarloResult{<:LineParameters}
    @test only(monte_carlo.details[:random].trials) == 3
    @test only(statistics(monte_carlo)) isa RLCG
    @test only(samples(monte_carlo)) isa RLCG
    @test only(histograms(monte_carlo)) isa RLCG
    @test first(only(histograms(monte_carlo)).R) isa HistogramDensity
    @test !applicable(Measurements.measurement, monte_carlo)

    line_frames=DataFrame(monte_carlo)
    @test size(line_frames) == size(only(statistics(monte_carlo)).R)
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

    parameterized=ParametricProblem(problem)
    @test_throws DomainError compute(parameterized, Combinatorial(Formulation()))
    skipped=compute(
        parameterized,
        Combinatorial(Formulation(); invalid = :skip)
    )
    @test skipped isa ParametricResult
    @test isempty(skipped)
    failures=skipped.details[:failures].values
    @test length(failures) == 2
    @test all(failure -> failure.exception_type == "DomainError", failures)
end

@testitem "DataModel / DataFrame / eager base-state presentation" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    design=only(Gridspace(PB.CableBuilder(
        "presentation-cable",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )))

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
    base_parameters=DataFrame(design, :baseparams)
    @test base_parameters.parameter == ["R", "L", "C"]
    @test base_parameters.unit == ["Ω/m", "H/m", "F/m"]
    @test_throws ErrorException DataFrame(design, :unsupported)
    constants=compute(CableConstantsProblem(design), Formulation())
    rendered=DataFrame(constants)
    @test rendered.parameter == ["R", "L", "C"]
    @test rendered.value == [constants.R, constants.L, constants.C]
end
