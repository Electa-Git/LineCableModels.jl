@testitem "ParametricBuilder and UQ / public API / calculation ownership" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine
    const ImportExport=LineCableModels.ImportExport
    const UQ=LineCableModels.UQ

    @test LineCableModels.AbstractProblemDefinition === Grammar.AbstractProblemDefinition
    @test LineCableModels.AbstractFormulation === Grammar.AbstractFormulation
    @test LineCableModels.AbstractProblemResult === Grammar.AbstractProblemResult
    @test LineCableModels.AbstractParametricResult === Grammar.AbstractParametricResult
    @test LineCableModels.AbstractUncertaintyResult === Grammar.AbstractUncertaintyResult
    @test LineCableModels.compute === Grammar.compute === Engine.compute
    @test LineCableModels.validate === LineCableModels.Validation.validate
    @test parentmodule(LineCableModels.validate) === LineCableModels.Validation
    @test length(methods(LineCableModels.validate)) == 1
    @test LineCableModels.observables === Grammar.observables
    @test LineCableModels.primitives === Grammar.primitives
    @test LineCableModels.preprocess === Grammar.preprocess
    @test parentmodule(Grammar.compute) === Grammar
    @test parentmodule(Grammar.observables) === Grammar
    @test parentmodule(Grammar.primitives) === Grammar
    @test parentmodule(Grammar.preprocess) === Grammar
    @test isempty(methods(Grammar.primitives))
    @test isempty(methods(Grammar.preprocess))
    @test LineCableModels.FormulationOptions === Grammar.FormulationOptions === NamedTuple
    @test LineCableModels.ComputationOptions === Grammar.ComputationOptions === NamedTuple
    @test LineCableModels.formulation_options === Grammar.formulation_options
    @test LineCableModels.computation_options === Grammar.computation_options
    for name in (
        :Grid, :AbsoluteError, :DeterministicGrid, :RelativeGrid, :AbsoluteGrid,
        :AbstractGrid, :AbstractUncertainGrid, :UncertainValue,
        :Gridspace
    )
        @test getproperty(LineCableModels, name) === getproperty(PB, name)
        @test parentmodule(getproperty(PB, name)) === PB
    end
    for name in (:Combinatorial, :ParametricProblem, :ParametricResult)
        @test getproperty(LineCableModels, name) === getproperty(PB, name)
        @test parentmodule(getproperty(PB, name)) === PB
    end
    for name in (
        :LinearError, :MonteCarlo, :LinearErrorResult, :MonteCarloResult,
        :SampleSummary, :HistogramDensity, :RLCG
    )
        @test getproperty(LineCableModels, name) === getproperty(UQ, name)
        @test parentmodule(getproperty(UQ, name)) === UQ
    end
    @test getproperty(LineCableModels, Symbol("@gridspace")) ===
          getproperty(PB, Symbol("@gridspace"))
    @test getproperty(LineCableModels, Symbol("@relax")) ===
          getproperty(PB, Symbol("@relax"))
    @test !isdefined(LineCableModels, :Computation)
    @test !isdefined(LineCableModels, Symbol("compute!"))
    @test !isdefined(LineCableModels, :CalculationManifest)
    @test !isdefined(LineCableModels, :ConfigurationFailure)
    @test !isdefined(LineCableModels, :manifest)
    @test !isdefined(PB, :CalculationManifest)
    @test !isdefined(PB, :ConfigurationFailure)
    @test !isdefined(PB, :manifest)
    @test !isdefined(Grammar, :Gridspace)
    @test !isdefined(Grammar, :MonteCarlo)
    @test Engine.AbstractProblemDefinition === Grammar.AbstractProblemDefinition
    @test Engine.AbstractFormulation === Grammar.AbstractFormulation
    @test parentmodule(Engine.AnalyticalInput) === Engine
    @test !isdefined(Engine, :EMTInput)
    @test !isdefined(ImportExport, :read_data)

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
    @test PB.PartBuilder !== nothing
    @test PB.PositionDefinition !== nothing
    @test !isdefined(PB, :MaterialDefinition)
    @test !isdefined(PB, :PartDefinition)
    @test !isdefined(PB, :CableDesignDefinition)
    @test !isdefined(PB, :EarthDefinition)
    @test !isdefined(PB, :SystemDefinition)
    @test !isdefined(PB, :_AbstractDefinition)
    @test isdefined(LineCableModels.DataModel, :Tubular)
    @test isdefined(LineCableModels.Materials, :Material)
    @test !isdefined(LineCableModels, :Configuration)
end

@testitem "Engine / grammar / strict analytical and computation options" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine
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

    @test !isdefined(Engine, :_FormulationFamily)
    @test !isdefined(Engine, :_ACTIVE_FORMULATION_BACKENDS)
    @test !isdefined(Engine, :_active_formulation_backend)
    @test !isdefined(Engine, :_activate_formulation_backend!)
    @test Formulation(:analytical) isa AnalyticalFormulation
    @test Formulation() isa AnalyticalFormulation

    options=LineCableModels.computation_options(Val(AnalyticalFormulation), (
        verbosity = (default = 1, NLsolve = 0),
        output_basis = :total
    ))
    @test options isa ComputationOptions
    @test options.output_basis == Val(:total)
    @test Engine.verbosity(options, :NLsolve) == 0
    @test Engine.verbosity(options, :unlisted) == 1
    @test_throws ArgumentError LineCableModels.computation_options(
        Val(AnalyticalFormulation),
        (unknown = true,)
    )
    @test_throws ArgumentError LineCableModels.computation_options(
        Val(AnalyticalFormulation),
        (output_basis = :unknown,)
    )
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
    @test [(item.layers[end].rho, item.layers[end].eps_r) for item in earth] == [
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
    system_space=systems
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

@testitem "ParametricBuilder and UQ / compute / deterministic and uncertainty ownership" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Random
    using Distributions
    import LineCableModels.ParametricBuilder as PB
    import LineCableModels.UQ as UQ

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    uncertain_radius=Grid((0.010, 0.012), 2.0)
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
    @test fieldnames(typeof(direct)) == (:formulation, :values)
    direct_observables=observables(direct)
    @test keys(direct_observables) == (:result,)
    @test direct_observables.result === result(direct)

    fixed_design=PB.CableBuilder(
        "fixed-compute-cable",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )
    fixed_problem=CableConstantsProblem(fixed_design)
    singleton_space=PB.Gridspace{CableConstantsProblem}(
        identity,
        (Grid(fixed_problem),)
    )
    singleton=compute(
        ParametricProblem(singleton_space),
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
    @test monte_carlo.trial_counts == [12, 12]
    @test length(samples(monte_carlo)) == 2
    @test length(histograms(monte_carlo)) == 2
    @test monte_carlo.root_seed == UInt64(0x1234)
    monte_carlo_observables=observables(monte_carlo)
    @test keys(monte_carlo_observables) ==
          (:result, :statistics, :samples, :histograms)
    @test monte_carlo_observables.result === result(monte_carlo)
    @test monte_carlo_observables.statistics === statistics(monte_carlo)
    @test monte_carlo_observables.samples === samples(monte_carlo)
    @test monte_carlo_observables.histograms === histograms(monte_carlo)

    resistance_pdf=first(histograms(monte_carlo)).R
    @test cdf(resistance_pdf, maximum(resistance_pdf)) == 1.0
    @test pdf(resistance_pdf, quantile(resistance_pdf, 0.5)) >= 0
    pdf_draw=rand(MersenneTwister(8), resistance_pdf)
    @test first(resistance_pdf.edges) <= pdf_draw <= last(resistance_pdf.edges)

    repeated=compute(ParametricProblem(space), policy)
    @test repeated.root_seed == monte_carlo.root_seed
    @test repeated.point_seeds == monte_carlo.point_seeds
    @test repeated.trial_counts == monte_carlo.trial_counts
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
    @test distribution_run.root_seed == UInt64(7)
    constants_frame=DataFrame(distribution_run)
    @test constants_frame.quantity == ["R", "L", "C"]
    @test constants_frame.trials == fill(4, 3)
    @test constants_frame.cdf_tol == fill(0.02, 3)
    @test !(:ci_half in propertynames(constants_frame))
    @test DataFrames.metadata(constants_frame, "monte_carlo").root_seed == UInt64(7)
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
    @test only(automatic_run.trial_counts) ==
          UQ._dkw_trials(3, 0.5, 0.9)

    propagated=compute(ParametricProblem(space), LinearError(formulation))
    @test propagated isa LinearErrorResult{<:CableConstants}
    @test fieldnames(typeof(propagated)) == (:formulation, :values)
    propagated_observables=observables(propagated)
    @test keys(propagated_observables) == (:result,)
    @test Measurements.cov(first(propagated).R, first(propagated).L) != 0
    @test !applicable(Measurements.measurement, monte_carlo)

    @test UQ._dkw_trials(3, 0.95, 0.02) ==
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

    deterministic_space=deterministic_problem
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

    uncertain_space=uncertain_problem
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
    @test only(monte_carlo.trial_counts) == 3
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

@testitem "ParametricBuilder / invalid points / natural failure" tags=[:unit] setup=[
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
    @test_throws MethodError Combinatorial(Formulation(); invalid = :skip)
    @test_throws MethodError LinearError(Formulation(); invalid = :skip)
    @test_throws MethodError MonteCarlo(Formulation(); invalid = :skip)
end

@testitem "DataModel / DataFrame / eager base-state presentation" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(; rho = 1.7241e-8)
    xlpe=PB.Material(; rho = 1.0e14, eps_r = 2.3)
    design=PB.CableBuilder(
        "presentation-cable",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(:core; thickness = 0.004, material = xlpe),
        PB.Conductor.Tubular(:screen; thickness = 0.001, material = copper),
        PB.Insulator.Tubular(:screen; thickness = 0.002, material = xlpe)
    )

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
