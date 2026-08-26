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
    @test LineCableModels.AbstractCoreResult === Grammar.AbstractCoreResult
    @test LineCableModels.AbstractResultSpace === Grammar.AbstractResultSpace
    @test LineCableModels.AbstractParametricResult === Grammar.AbstractParametricResult
    @test LineCableModels.AbstractUncertaintyResult === Grammar.AbstractUncertaintyResult
    @test LineCableModels.compute === Grammar.compute === Engine.compute
    @test LineCableModels.validate === LineCableModels.Validation.validate
    @test parentmodule(LineCableModels.validate) === LineCableModels.Validation
    @test length(methods(LineCableModels.validate)) == 1
    @test LineCableModels.observables === Grammar.observables
    @test LineCableModels.observe === Grammar.observe
    @test LineCableModels.computation_details === Grammar.computation_details
    @test LineCableModels.details === Grammar.details
    @test parentmodule(Grammar.compute) === Grammar
    @test parentmodule(Grammar.computation_details) === Grammar
    @test parentmodule(Grammar.details) === Grammar
    @test parentmodule(Grammar.observables) === Grammar
    @test parentmodule(Grammar.observe) === Grammar
    @test parentmodule(Grammar.check_core_result) === Grammar
    @test Base.ispublic(Grammar, :check_core_result)
    @test !isdefined(LineCableModels, :check_core_result)
    @test LineCableModels.FormulationOptions === Grammar.FormulationOptions === NamedTuple
    @test LineCableModels.ComputationOptions === Grammar.ComputationOptions === NamedTuple
    @test LineCableModels.ComputationDetails === Grammar.ComputationDetails === NamedTuple
    @test LineCableModels.formulation_options === Grammar.formulation_options
    @test LineCableModels.computation_options === Grammar.computation_options
    @test getproperty(LineCableModels, Symbol("@observe")) ===
          getproperty(Grammar, Symbol("@observe"))
    struct UnregisteredDetailsOwner end
    @test_throws MethodError Grammar.computation_details(
        Val(UnregisteredDetailsOwner),
        nothing
    )
    @test_throws MethodError Grammar.computation_details(nothing)
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
        :SampleSummary, :HistogramDensity
    )
        @test getproperty(LineCableModels, name) === getproperty(UQ, name)
        @test parentmodule(getproperty(UQ, name)) === UQ
    end
    for name in (
        :root_seed, :point_seed, :trial_count,
        :confidence, :cdf_tolerance, :sampling_distribution
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
    @test !isdefined(LineCableModels, :primitives)
    @test !isdefined(LineCableModels, :preprocess)
    @test !isdefined(PB, :CalculationManifest)
    @test !isdefined(PB, :ConfigurationFailure)
    @test !isdefined(PB, :manifest)
    @test !isdefined(Grammar, :Gridspace)
    @test !isdefined(Grammar, :MonteCarlo)
    @test !isdefined(Grammar, :primitives)
    @test !isdefined(Grammar, :preprocess)
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
    @test !isdefined(LineCableModels, :RLC)
    @test !isdefined(LineCableModels, :RLCG)
    @test !isdefined(UQ, :RLC)
    @test !isdefined(UQ, :RLCG)
end

@testitem "Grammar / results / core and result-space taxonomy" tags=[:unit] begin
    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ
    const Engine=LineCableModels.Engine

    abstract type AbstractExternalPayload end
    struct ExternalPayload <: AbstractExternalPayload
        value::Int
    end

    @test LineCableModels.CableConstants <: Grammar.AbstractCoreResult
    @test LineCableModels.LineParameters <: Grammar.AbstractCoreResult
    @test !(Engine.LineParametersTrace <: Grammar.AbstractCoreResult)
    @test !(Engine.LineParametersBenchmark <: Grammar.AbstractCoreResult)
    @test PB.ParametricResult <: Grammar.AbstractResultSpace
    @test UQ.LinearErrorResult <: Grammar.AbstractResultSpace
    @test UQ.MonteCarloResult <: Grammar.AbstractResultSpace

    @test Grammar.check_core_result(ExternalPayload) === nothing
    external=PB.ParametricResult(nothing, ExternalPayload[ExternalPayload(1)])
    @test only(external).value == 1

    nested=PB.ParametricResult(nothing, ExternalPayload[ExternalPayload(2)])
    @test_throws ArgumentError PB.ParametricResult(nothing, [nested])
    @test_throws ArgumentError PB.ParametricResult(
        nothing,
        AbstractExternalPayload[ExternalPayload(3)]
    )
    @test_throws ArgumentError PB.ParametricResult(nothing, Any[ExternalPayload(4)])
end

@testitem "Grammar / results / standard collection semantics" tags=[:unit] begin
    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct CollectionPayload
        value::Int
    end

    expected=CollectionPayload[CollectionPayload(1), CollectionPayload(2)]
    spaces=(
        PB.ParametricResult(nothing, copy(expected)),
        UQ.LinearErrorResult(nothing, copy(expected)),
        UQ.MonteCarloResult(
            nothing,
            copy(expected),
            [:first, :second],
            nothing,
            nothing,
            UInt64(7),
            UInt64[8, 9],
            Int[10, 10]
        )
    )

    for space in spaces
        @test Base.IteratorSize(typeof(space)) == Base.HasShape{1}()
        @test Base.IteratorEltype(typeof(space)) == Base.HasEltype()
        @test eltype(space) === CollectionPayload
        @test length(space) == 2
        @test size(space) == (2,)
        @test firstindex(space) == 1
        @test lastindex(space) == 2
        @test space[2] === expected[2]
        @test first(space) === expected[1]
        @test last(space) === expected[2]
        @test collect(space) == expected
        @test map(value -> value.value, space) == [1, 2]
        @test collect(zip(space, (:a, :b))) == [(expected[1], :a), (expected[2], :b)]
        @test_throws ArgumentError only(space)
    end

    singleton=PB.ParametricResult(nothing, CollectionPayload[expected[1]])
    @test only(singleton) === first(singleton)
    @test all(
        method -> method.module ∉ (LineCableModels, Grammar, PB, UQ),
        methods(Base.only)
    )
end

@testitem "ParametricBuilder and UQ / compute / empty problem spaces" tags=[:unit] begin
    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct EmptyFormulation <: Grammar.AbstractFormulation end

    space=PB.Gridspace{Int}(identity, (PB.Grid(()),))
    problem=PB.ParametricProblem(space)
    @test isempty(space)
    @test_throws ArgumentError Grammar.compute(
        problem,
        PB.Combinatorial(EmptyFormulation())
    )
    @test_throws ArgumentError Grammar.compute(problem, UQ.LinearError(EmptyFormulation()))
    @test_throws ArgumentError Grammar.compute(
        problem,
        UQ.MonteCarlo(EmptyFormulation(); trials = 2, seed = 1)
    )
end

@testitem "ParametricBuilder and UQ / compute / exact traversal grammars" tags=[:unit] begin
    using Statistics

    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct TraversalProblem <: Grammar.AbstractProblemDefinition
        value::Int
    end
    struct StableTraversalFormulation <: Grammar.AbstractFormulation end
    struct UnstableTraversalFormulation <: Grammar.AbstractFormulation end
    struct DetailsTraversalFormulation <: Grammar.AbstractFormulation end
    struct TraversalResult
        value::Int
    end

    Grammar.compute(problem::TraversalProblem, ::StableTraversalFormulation;
        options = (;)) = TraversalResult(problem.value)
    function Grammar.compute(
            problem::TraversalProblem,
            ::UnstableTraversalFormulation;
            options = (;)
    )
        return problem.value == 1 ? TraversalResult(problem.value) : problem.value
    end
    Grammar.compute(problem::TraversalProblem, ::DetailsTraversalFormulation;
        options = (;)) = TraversalResult(problem.value)
    function Grammar.computation_details(
            ::Val{DetailsTraversalFormulation},
            output::TraversalResult
    )::Grammar.ComputationDetails
        return output.value == 1 ? (raw = Dict("value" => 1),) : (raw = [2],)
    end

    space=PB.Gridspace{TraversalProblem}(
        TraversalProblem,
        (PB.Grid((1, 2)),)
    )
    problem=PB.ParametricProblem(space)
    stable=Grammar.compute(
        problem,
        PB.Combinatorial(StableTraversalFormulation())
    )
    @test stable.values isa Vector{TraversalResult}
    @test getproperty.(stable, :value) == [1, 2]
    @test_throws ArgumentError Grammar.compute(
        problem,
        PB.Combinatorial(UnstableTraversalFormulation())
    )
    @test_throws ArgumentError Grammar.compute(
        problem,
        PB.Combinatorial(
            DetailsTraversalFormulation();
            options = (retain_details = true,)
        )
    )

    struct MonteTraversalFormulation <: Grammar.AbstractFormulation end
    struct MonteTraversalResult
        value::Float64
    end
    Grammar.compute(problem::TraversalProblem, ::MonteTraversalFormulation;
        options = (;)) = MonteTraversalResult(problem.value)
    UQ._observable_count(::MonteTraversalResult) = 1
    UQ._sample_storage(::MonteTraversalResult, trials::Int) = Vector{Float64}(undef, trials)
    UQ._sample_axis(::MonteTraversalResult) = nothing
    function UQ._record_sample!(
            storage::Vector{Float64},
            value::MonteTraversalResult,
            trial::Int,
            ::Nothing
    )
        storage[trial] = value.value
        return storage
    end
    function UQ._aggregate(
            storage::Vector{Float64},
            first_result::MonteTraversalResult,
            ::UQ.MonteCarlo
    )
        summary = UQ.SampleSummary(storage)
        product = first_result.value == 1 ? summary : (summary = summary,)
        return (
            representation = MonteTraversalResult(mean(storage)),
            statistics = product,
            samples = nothing,
            histograms = nothing
        )
    end
    @test_throws ArgumentError Grammar.compute(
        problem,
        UQ.MonteCarlo(MonteTraversalFormulation(); trials = 2, seed = 3)
    )

    @test Grammar.computation_owner(StableTraversalFormulation()) ===
          StableTraversalFormulation
    @test Base.ispublic(Grammar, :computation_owner)
    @test Base.ispublic(Grammar, :detach)
    @test Base.ispublic(PB, :traverse)
    @test Base.ispublic(PB, :sample_uncertainty)
    @test Base.ispublic(LineCableModels.ReportBuilder, :clip)
    @test Base.ispublic(LineCableModels.ReportBuilder, :encode_cell)
    @test Base.ispublic(LineCableModels.ReportBuilder, :XLSXSheet)
    @test Base.ispublic(LineCableModels.ReportBuilder, :XLSXWorkbook)
    @test Base.ispublic(LineCableModels.PlotBuilder, :validate_export_theme)
    for name in (
        :computation_owner, :traverse, :sample_uncertainty,
        :clip, :encode_cell, :validate_export_theme
    )
        @test !isdefined(LineCableModels, name)
    end
    @test LineCableModels.detach === Base.detach
    @test LineCableModels.detach !== Grammar.detach
    @test !Base.isexported(LineCableModels, :detach)

    source_root=joinpath(pkgdir(LineCableModels), "src", "parametricbuilder")
    @test isfile(joinpath(source_root, "traversal.jl"))
    @test isfile(joinpath(source_root, "engine", "cableconstants.jl"))
    @test !isfile(joinpath(source_root, "compute.jl"))
    for name in (:_traverse, :_computation_owner, :_append_result, :_sample_uncertainty)
        @test !isdefined(PB, name)
    end
    @test !isdefined(Grammar, :_detach_and_scale)
    @test !isdefined(LineCableModels.ReportBuilder, :_clip_field)
    for private_name in (
        :_xlsx_string,
        :_xlsx_strings,
        :_xlsx_units,
        :_xlsx_destination,
        :_write_xlsx_sheet!
    )
        @test !isdefined(LineCableModels.ReportBuilder, private_name)
    end
    @test !isdefined(LineCableModels.PlotBuilder, :_validate_export_theme)
    @test !isdefined(LineCableModels.Engine, :_has_uncertainty_type)
    @test !isdefined(LineCableModels.Engine, :_cable_constants_problem)
    @test !isdefined(LineCableModels.DataModel, :_compute_cable_constants)
    @test !isdefined(LineCableModels.DataModel, :_base_parameters)
    @test !isdefined(LineCableModels.ImportExport, :_serialize_value)
    @test !isdefined(LineCableModels.ImportExport, :_deserialize_value)
    @test !isdefined(LineCableModels.ImportExport, :_deserialize_extension)
end

@testitem "ParametricBuilder and UQ / supplemental output / external owner retention" tags=[:unit] begin
    using Statistics

    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct SupplementalProblem <: Grammar.AbstractProblemDefinition
        value::Float64
    end
    struct SupplementalFormulation <: Grammar.AbstractFormulation end
    struct SupplementalResult <: Grammar.AbstractProblemResult
        value::Float64
        raw::Dict{String, Any}
    end
    const supplemental_detail_calls=Ref(0)

    function Grammar.compute(
            problem::SupplementalProblem,
            ::SupplementalFormulation;
            options::NamedTuple = (;)
    )
        isempty(options) || throw(ArgumentError("supplemental test options must be empty"))
        return SupplementalResult(
            problem.value,
            Dict{String, Any}("coordinate" => problem.value)
        )
    end

    function Grammar.computation_details(
            ::Val{SupplementalFormulation},
            output::SupplementalResult
    )::Grammar.ComputationDetails
        supplemental_detail_calls[]+=1
        return (
            diagnostics = (iterations = 1,),
            raw = output.raw
        )
    end

    UQ._observable_count(::SupplementalResult) = 1
    UQ._sample_storage(::SupplementalResult, trials::Int) = Vector{Float64}(undef, trials)
    UQ._sample_axis(::SupplementalResult) = nothing
    function UQ._record_sample!(
            storage::Vector{Float64},
            value::SupplementalResult,
            trial::Int,
            ::Nothing
    )
        storage[trial] = value.value
        return storage
    end
    function UQ._aggregate(
            storage::Vector{Float64},
            ::SupplementalResult,
            ::UQ.MonteCarlo
    )
        summary = UQ.SampleSummary(storage)
        representation = SupplementalResult(
            summary.mean,
            Dict{String, Any}("aggregate" => summary.mean)
        )
        return (
            representation,
            statistics = summary,
            samples = nothing,
            histograms = nothing
        )
    end

    formulation=SupplementalFormulation()
    deterministic_space=PB.Gridspace{SupplementalProblem}(
        SupplementalProblem,
        (PB.Grid((1.0, 2.0)),)
    )
    problem=PB.ParametricProblem(deterministic_space)

    ordinary=Grammar.compute(problem, PB.Combinatorial(formulation))
    @test Grammar.details(ordinary) === (;)
    @test fieldtype(typeof(ordinary), :details) === @NamedTuple{}
    @test supplemental_detail_calls[] == 0

    retained=Grammar.compute(
        problem,
        PB.Combinatorial(
            formulation;
            options = (retain_details = true,)
        )
    )
    @test @inferred(Grammar.details(retained)) === retained.details
    @test length(retained.details.points) == 2
    @test retained.details.points[1].diagnostics == (iterations = 1,)
    @test retained.details.points[2].raw["coordinate"] == 2.0
    @test fieldtype(typeof(retained), :details) <: NamedTuple
    @test supplemental_detail_calls[] == 2

    uncertain_space=PB.Gridspace{SupplementalProblem}(
        SupplementalProblem,
        (PB.Grid(1.0, 5.0),)
    )
    monte_problem=PB.ParametricProblem(uncertain_space)
    supplemental_detail_calls[]=0
    default_monte=Grammar.compute(
        monte_problem,
        UQ.MonteCarlo(formulation; trials = 3, seed = 41)
    )
    @test Grammar.details(default_monte) === (;)
    Grammar.details(default_monte)
    @test @allocated(Grammar.details(default_monte)) == 0
    @test supplemental_detail_calls[] == 0

    retained_monte=Grammar.compute(
        monte_problem,
        UQ.MonteCarlo(
            formulation;
            trials = 3,
            seed = 41,
            options = (retain_details = true,)
        )
    )
    @test length(retained_monte.details.trials) == 1
    @test length(only(retained_monte.details.trials)) == 3
    @test all(
        record -> keys(record) == (:diagnostics, :raw),
        only(retained_monte.details.trials)
    )
    @test retained_monte.root_seed == default_monte.root_seed
    @test retained_monte.point_seeds == default_monte.point_seeds
    @test retained_monte.trial_counts == default_monte.trial_counts
    @test only(result(retained_monte)).value == only(result(default_monte)).value
    @test supplemental_detail_calls[] == 3

    @test Grammar.computation_options(Val(PB.Combinatorial), (;)) ==
          (retain_details = false,)
    @test Grammar.computation_options(Val(UQ.LinearError), (;)) ==
          (retain_details = false,)
    @test Grammar.computation_options(Val(UQ.MonteCarlo), (;)) ==
          (retain_details = false,)
    @test_throws ArgumentError Grammar.computation_options(
        Val(PB.Combinatorial),
        (retain_details = :yes,)
    )
    @test_throws ArgumentError Grammar.computation_options(
        Val(UQ.MonteCarlo),
        (unknown = true,)
    )
    @test_throws MethodError Grammar.computation_options(
        Val(PB.Combinatorial),
        Dict(:retain_details => true)
    )
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

    options=LineCableModels.computation_options(Val(AnalyticalFormulation),
        (
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
    @test fieldnames(typeof(direct)) == (:formulation, :values, :details)
    @test @inferred(details(direct)) === (;)
    @test_throws MethodError observables(direct)

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
    retained_singleton=compute(
        ParametricProblem(singleton_space),
        Combinatorial(
            formulation;
            options = (retain_details = true,)
        )
    )
    @test retained_singleton.details == (points = [(;)],)

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
    @test all(product -> keys(product) == (:R, :L, :C), statistics(monte_carlo))
    @test all(product -> keys(product) == (:R, :L, :C), samples(monte_carlo))
    @test all(product -> keys(product) == (:R, :L, :C), histograms(monte_carlo))
    @test eltype(statistics(monte_carlo)) === typeof(first(statistics(monte_carlo)))
    @test eltype(samples(monte_carlo)) === typeof(first(samples(monte_carlo)))
    @test eltype(histograms(monte_carlo)) === typeof(first(histograms(monte_carlo)))
    @test !(first(statistics(monte_carlo)) isa CableConstants)
    @test_throws MethodError CableConstants([1.0], [2.0], [3.0])
    @test minimum(samples(monte_carlo)[1].R) >
          maximum(samples(monte_carlo)[2].R)
    @test monte_carlo.root_seed == UInt64(0x1234)
    @test @inferred(details(monte_carlo)) === (;)
    @test_throws MethodError observables(monte_carlo)
    @test observe(monte_carlo, statistics, R, mean, 1) ==
          mean(first(samples(monte_carlo)).R)
    @test observe(monte_carlo, samples, R, 1, :) ==
          first(samples(monte_carlo)).R
    @test observe(monte_carlo, histograms, R, 1) ===
          first(histograms(monte_carlo)).R
    @test !applicable(observe, first(statistics(monte_carlo)), R)
    multi_point_table=DataFrame(monte_carlo)
    @test multi_point_table isa DataFrame
    @test names(multi_point_table) == [
        "point", "quantity", "mean", "std", "min", "q05", "median", "q95",
        "max", "n", "unit", "trials", "confidence", "cdf_tol"
    ]
    @test multi_point_table.point == repeat(1:2; inner = 3)
    @test multi_point_table.quantity == repeat([:R, :L, :C], 2)
    monte_carlo_metadata=DataFrames.metadata(
        multi_point_table,
        "monte_carlo"
    )
    @test monte_carlo_metadata.root_seed == monte_carlo.root_seed
    @test monte_carlo_metadata.point_seeds == monte_carlo.point_seeds
    @test monte_carlo_metadata.trial_counts == monte_carlo.trial_counts
    @test monte_carlo_metadata.confidence == confidence(monte_carlo)
    @test monte_carlo_metadata.cdf_tol == cdf_tolerance(monte_carlo)
    @test monte_carlo_metadata.distribution == sampling_distribution(monte_carlo)
    @test DataFrames.metadata(multi_point_table, "units") == Dict(
        :R => "Ω/km",
        :L => "mH/km",
        :C => "μF/km"
    )

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
    @test constants_frame.quantity == [:R, :L, :C]
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
    @test fieldnames(typeof(propagated)) == (:formulation, :values, :details)
    @test @inferred(details(propagated)) === (;)
    @test_throws MethodError observables(propagated)
    @test Measurements.cov(first(propagated).R, first(propagated).L) != 0
    @test !applicable(Measurements.measurement, monte_carlo)

    @test UQ._dkw_trials(3, 0.95, 0.02) ==
          ceil(Int, log(2 * 3 / 0.05) / (2 * 0.02^2))
end

@testitem "Engine / Formulation / materialized and Gridspace line problems" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    import LineCableModels.ParametricBuilder as PB
    using Statistics

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
    @test keys(only(statistics(monte_carlo))) == (:R, :L, :C, :G)
    @test keys(only(samples(monte_carlo))) == (:R, :L, :C, :G)
    @test keys(only(histograms(monte_carlo))) == (:R, :L, :C, :G)
    @test first(only(histograms(monte_carlo)).R) isa HistogramDensity
    @test basis(monte_carlo) === :pul
    @test observe(monte_carlo, statistics, R, Statistics.mean, 1, :, :, :) ==
          Statistics.mean.(only(statistics(monte_carlo)).R)
    @test observe(monte_carlo, samples, R, 1, 1, 1, 1, :) ==
          only(samples(monte_carlo)).R[1, 1, 1, :]
    @test observe(monte_carlo, histograms, R, 1, 1, 1, 1) ===
          only(histograms(monte_carlo)).R[1, 1, 1]
    @test !applicable(Measurements.measurement, monte_carlo)

    line_table=DataFrame(monte_carlo)
    @test line_table isa DataFrame
    @test names(line_table) == [
        "point", "row", "column", "frequency", "quantity", "mean", "std",
        "min", "q05", "median", "q95", "max", "n", "unit", "trials",
        "confidence", "cdf_tol"
    ]
    @test DataFrames.nrow(line_table) ==
          4 * prod(size(only(statistics(monte_carlo)).R))
    @test unique(line_table.quantity) == [:R, :L, :C, :G]
    @test all(==("Ω/km"), line_table[line_table.quantity .== :R, :unit])
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
