@testitem "Gauntlet / one formulation provenance writer retains route inputs" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    ordinary = Formulation(insulation_admittance=:Ametani2004)
    record = formulation_record(ordinary)
    @test record.schema_version == 1
    @test record.insulation_admittance.identifier === :Ametani2004
    @test record.insulation_admittance.routes !== nothing
    @test record.equivalent_earth.rule === :Layer
    @test record.equivalent_earth.layer == -1
    @test basename(String(which(formulation_record, (typeof(ordinary),)).file)) == "provenance.jl"
    make_route(scale) = (material, frequency, temperature, assumptions) ->
        scale * constitutive(ordinary.methods.insulation_admittance, material, frequency, temperature)
    first_formulation = Formulation(insulation_admittance=formula(:Ametani2004; route=make_route(1.0)))
    second_formulation = Formulation(insulation_admittance=formula(:Ametani2004; route=make_route(2.0)))
    @test typeof(first_formulation.methods.insulation_admittance.route) ===
        typeof(second_formulation.methods.insulation_admittance.route)
    first_record = formulation_record(first_formulation)
    second_record = formulation_record(second_formulation)
    @test GauntletSupport.semantic_sha256(first_record) !=
        GauntletSupport.semantic_sha256(second_record)
    @test GauntletSupport._semantic_formulation_record(first_record) !=
        GauntletSupport._semantic_formulation_record(second_record)
    @test_throws ArgumentError GauntletSupport._semantic_formulation_record(
        (type="legacy.LineParametersFormulation",))
    @test_throws ArgumentError GauntletSupport._selection_value(Ref(1.0))
    @test GauntletSupport._same_problem_structure((owner=Group,), (owner=Group,))
    @test !GauntletSupport._same_problem_structure((owner=Group,), (owner=Assembly,))
    @test GauntletSupport._same_problem_structure(Disk{Float64}, Disk{Float64})
    @test !GauntletSupport._same_problem_structure(Disk{Float64}, Disk{Float32})
end

@testitem "Gauntlet / resolved geometry and flattening participate in reuse" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor, rho=1.72e-8)
    design = build(CableDesign, "resolved-provenance",
        terminal(:core, solid(copper, Disk(0.005))))
    system = build(LineCableSystem, design, Pose2(0.0, -0.1);
        connections=Dict(:core=>1))
    problem = LineParametersProblem(system; frequencies=[50.0],
        earth_props=homogeneous(rho=100.0))
    declaration = LineCableModels.ImportExport.serialize_value(problem)
    original = numerical_input_sha256(problem)
    @test original == numerical_input_sha256(deepcopy(problem))
    restored = LineCableModels.ImportExport.deserialize_value(declaration)
    @test original == numerical_input_sha256(restored)
    # A later resolver may change a physical region without changing the
    # declaration. Do not amend the JSON transport schema to fingerprint it.
    placed = first(design.geometry.regions)
    design.geometry.regions[1] = DM.PlacedRegion(placed.source,
        Disk(placed.primitive.r * 0.9, placed.primitive.at), placed.terminal,
        placed.placement, placed.paths)
    @test LineCableModels.ImportExport.serialize_value(problem) == declaration
    @test numerical_input_sha256(problem) != original
    paths = getproperty.(implementation_record(Formulation()).blobs, :path)
    @test all(path -> path in paths, (
        "src/datamodel/baseparams/geometry.jl",
        "src/datamodel/baseparams/resistance.jl",
        "src/datamodel/baseparams/inductance.jl",
        "src/datamodel/baseparams/dielectrics.jl",
        "src/datamodel/placement/paths.jl"))
    @test !any(path -> startswith(path, "src/plotbuilder/") ||
        startswith(path, "ext/LineCableModelsMakieExt/"), paths)
    @test !("src/engine/earthimpedance/formulas/carson1926.jl" in paths)
end
