@testitem "Gauntlet / comparison entry point resolves provenance and live catalogue" tags=[:gauntlet_toolkit] setup=[TestFixtures] begin
    using .TestFixtures
    include(joinpath(pkgdir(LineCableModels), "test", "gauntlet", "formulation_comparisons.jl"))
    records=catalogue()
    expected=Set(vcat(
        [(:earth_impedance, id) for id in EARTH_IMPEDANCE.formulas()],
        [(:earth_admittance, id) for id in EARTH_ADMITTANCE.formulas()]))
    @test Set((record.kind, record.identifier) for record in records) == expected
    selected=variant(first(records))
    model=(nominal_problem = TestFixtures.line_parameters_problem(; frequencies = [50.0]),)
    provenance=comparison_provenance(model, selected)
    @test !isempty(provenance.implementation.blobs)
    @test provenance.input_sha256 == numerical_input_sha256(model.nominal_problem)
    @test provenance.implementation.selection.schema_version == 2
end

@testitem "Gauntlet / PSCAD reference entry point uses native catalogue" tags=[:gauntlet_toolkit] begin
    # A standalone entry point owns its backend module. Loading it into the
    # shared test process would replace the public :pscad constructor methods.
    script = raw"""
    using Test, LineCableModels
    include(joinpath(pkgdir(LineCableModels), "test", "gauntlet", "pscad_reference.jl"))
    const P = GauntletSupport.PSCADBenchmarks
    @testset "PSCAD reference catalogue" begin
        expected = Set((placement, identifier)
            for placement in (:overhead, :underground, :mixed)
            for identifier in P.formulas(Val(placement)))
        @test Set((record.field, record.id) for record in PSCAD_CATALOGUE) == expected
        @test allunique(variant_id.(PSCAD_CATALOGUE))
        for record in PSCAD_CATALOGUE
            @test record.selector === record.id
            @test P.pscad_setting(Val(record.selector), Val(record.field)).value >= 0
        end
    end
    """
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`
    @test success(command)
end

@testitem "Gauntlet / catalogue applicability follows formula-owned preflight" tags=[:gauntlet_toolkit] setup=[TestFixtures] begin
    using .TestFixtures
    include(joinpath(pkgdir(LineCableModels), "test", "gauntlet", "fem_catalogue.jl"))
    const EP=LineCableModels.Earth
    problem=TestFixtures.line_parameters_problem(; frequencies = [50.0])
    model=(; problem)
    workspace=prepare_case(model)
    choose(kind, id)=variant(coverage_record(kind, id))
    reason(kind, id)=case_skip_reason(model, choose(kind, id), workspace)
    @test isnothing(reason(:earth_impedance, :default))
    @test occursin("formula not implemented for source in layer", reason(:earth_impedance, :Ametani2009))
    @test occursin("formula not implemented for source in layer", reason(:earth_impedance, :Carson1926))

    layered=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 4.0), EP.EarthLayer(200.0, 20.0, 1.0)))
    layered_model=(problem = LineParametersProblem(problem.system;
        earth_props = layered, frequencies = [50.0]),)
    layered_input=prepare_case(layered_model)
    @test !isnothing(case_skip_reason(layered_model, choose(:earth_impedance, :default), layered_input))

    overhead_system=build(LineCableSystem, problem.system.designs,
        [Pose2(p.x, -p.y) for p in problem.system.positions];
        connections = problem.system.connections, system_id = "overhead-preflight")
    overhead_model=(problem = LineParametersProblem(overhead_system;
        earth_props = problem.earth_props, frequencies = [50.0]),)
    overhead_workspace=prepare_case(overhead_model)
    @test isnothing(case_skip_reason(overhead_model,
        choose(:earth_impedance, :default), overhead_workspace))
    @test isnothing(case_skip_reason(overhead_model,
        choose(:earth_impedance, :Carson1926), overhead_workspace))
    @test occursin("formula not implemented for source in layer",
        case_skip_reason(overhead_model,
            choose(:earth_impedance, :Saad1996), overhead_workspace))

    vertical_system=build(LineCableSystem, problem.system.designs,
        [Pose2(0.0, -Float64(i)) for i in eachindex(problem.system.designs)];
        connections = problem.system.connections, system_id = "vertical-pair-preflight")
    vertical_model=(problem = LineParametersProblem(vertical_system;
        earth_props = problem.earth_props, frequencies = [50.0]),)
    vertical_workspace=prepare_case(vertical_model)
    @test occursin("nonzero horizontal",
        case_skip_reason(vertical_model,
            choose(:earth_impedance, :Saad1996), vertical_workspace))
    @test isnothing(case_skip_reason(vertical_model,
        choose(:earth_impedance, :default), vertical_workspace))

    calls=Ref(0)
    failing_kernel=(
        functor, pair, workspace)->(calls[]+=1; throw(DomainError(0, "numerical sentinel")))
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{
            :default, typeof(LineCableModels.Engine.EarthImpedance.earth_impedance)},
        ::$(typeof(failing_kernel))) = (;)
    overridden=merge(choose(:earth_impedance, :default),
        (earth_impedance = formula(:default; hooks = (contribution = failing_kernel,)),))
    @test isnothing(case_skip_reason(model, overridden, workspace))
    @test calls[] == 0
    @test_throws DomainError compute(problem, formulation(overridden))
    @test calls[] == 1

    # A single design can contain several independent coaxial assemblies.
    copper=TestFixtures.copper_material()
    insulator=TestFixtures.insulator_material()
    member=terminal(:core, solid(copper, Disk(1e-3)), insulation(insulator; t = 1e-3))
    design=build(CableDesign, "three-assemblies", cores(member;
        n = 3, r = 0.01, names = (:a, :b, :c)))
    system=build(LineCableSystem, [design], [Pose2(0.0, -1.0)];
        connections = [Dict(:a=>1, :b=>2, :c=>3)], system_id = "single-design-many-assemblies")
    multiple=(problem = LineParametersProblem(system;
        earth_props = problem.earth_props, frequencies = [50.0]),)
    many_workspace=prepare_case(multiple)
    @test length(system.designs) == 1
    @test length(many_workspace.cable.assemblies) == 3
    @test occursin("formula not implemented for source in layer",
        case_skip_reason(multiple,
            choose(:earth_impedance, :Ametani2009), many_workspace))

    # An author-independent catalogue description must not claim case applicability.
    for record in catalogue()
        @test !hasproperty(record, :applicable)
        @test !hasproperty(record, :category)
    end
end

@testitem "Gauntlet / manual PSCAD benchmark selections use registered formulations" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport
    const P=GauntletSupport.PSCADBenchmarks
    directory=joinpath(pkgdir(LineCableModels), "test", "gauntlet", "benchmarks", "pscad")
    paths=filter(path->startswith(basename(path), "benchmark_")&&endswith(path, ".jl"),
        readdir(directory; join = true))
    @test !isempty(paths)
    for path in sort(paths)
        @testset "$(basename(path))" begin
            parsed = Meta.parseall(read(path, String))
            blocks = filter(
                node -> node isa Expr && node.head === :macrocall &&
                        first(node.args) === Symbol("@testitem"),
                parsed.args)
            @test !isempty(blocks)
            for block in blocks
                declarations = Dict(node.args[1] => node.args[2]
                for node in last(block.args).args
                if node isa Expr && node.head === :(=) &&
                   node.args[1] in (:reference_formulation, :candidate_formulation))
                for (name, expected) in ((:reference_formulation, P.PSCADFormulation),
                    (:candidate_formulation, LineParametersFormulation))
                    @test haskey(declarations, name)
                    expression = declarations[name]
                    # Evaluate only the existing native construction calls.
                    # Never include the manual body: it can run external solvers.
                    constructor_call = expression isa Expr && expression.head === :call &&
                                       first(expression.args) === :Formulation
                    @test constructor_call
                    constructor_call || continue
                    value = Core.eval(@__MODULE__, expression)
                    for selected in (value isa Gridspace ? value : (value,))
                        @test selected isa expected
                    end
                end
            end
        end
    end
end
