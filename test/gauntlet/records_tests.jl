@testitem "Gauntlet / one formulation record writer retains route inputs" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport.Gauntlet
    ordinary=Formulation(insulation_admittance = :Ametani2004)
    record=formulation_record(ordinary)
    @test record.schema_version == 2
    @test record.insulation_admittance.identifier === :Ametani2004
    @test record.insulation_admittance.binding !== nothing
    @test record.earth_impedance.equivalent_earth === nothing
    @test record.earth_admittance.equivalent_earth === nothing
    @test basename(String(which(formulation_record, (typeof(ordinary),)).file)) ==
          "records.jl"
    make_route(scale)=(material,
        frequency,
        temperature,
        assumptions, options,
        workspace)->scale*constitutive(
        ordinary.methods.insulation_admittance, material, frequency, temperature)
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{:Ametani2004,
            typeof(LineCableModels.Engine.InsulationAdmittance.insulation_material)},
        ::$(typeof(make_route(1.0)))) = (;)
    first_formulation=Formulation(insulation_admittance = formula(
        :Ametani2004; hooks = (contribution = make_route(1.0),)))
    second_formulation=Formulation(insulation_admittance = formula(
        :Ametani2004; hooks = (contribution = make_route(2.0),)))
    @test typeof(first_formulation.methods.insulation_admittance.hooks.contribution) ===
          typeof(second_formulation.methods.insulation_admittance.hooks.contribution)
    first_record=formulation_record(first_formulation)
    second_record=formulation_record(second_formulation)
    @test GauntletSupport.Gauntlet.semantic_sha256(first_record) !=
          GauntletSupport.Gauntlet.semantic_sha256(second_record)
    @test_throws ArgumentError GauntletSupport.Gauntlet._selection_value(Ref(1.0))
    @test GauntletSupport.Gauntlet._selection_value(Disk{Float64}) ==
          sprint(show, Disk{Float64}; context = (:module=>nothing, :compact=>false))
    @test GauntletSupport.Gauntlet._selection_value(Disk(1.0)).type ==
          GauntletSupport.Gauntlet._selection_value(typeof(Disk(1.0)))
end

@testitem "Gauntlet / resolved geometry and flattening participate in reuse" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport.Gauntlet
    const DM=LineCableModels.DataModel
    copper=Material(kind = :conductor, rho = 1.72e-8)
    design=build(CableDesign, "resolved-geometry",
        terminal(:core, solid(copper, Disk(0.005))))
    system=build(LineCableSystem, design, Pose2(0.0, -0.1);
        connections = Dict(:core=>1))
    problem=LineParametersProblem(system; frequencies = [50.0],
        earth_props = homogeneous(rho = 100.0))
    declaration=LineCableModels.ImportExport.serialize_value(problem)
    original=numerical_input_sha256(problem)
    @test original == numerical_input_sha256(deepcopy(problem))
    restored=LineCableModels.ImportExport.deserialize_value(declaration)
    @test original == numerical_input_sha256(restored)
    # A later resolver may change a physical region without changing the
    # declaration. Do not amend the JSON transport schema to fingerprint it.
    placed=first(design.geometry.regions)
    design.geometry.regions[1]=DM.PlacedRegion(placed.source,
        Disk(placed.primitive.r*0.9, placed.primitive.at), placed.terminal,
        placed.placement, placed.paths)
    @test LineCableModels.ImportExport.serialize_value(problem) == declaration
    @test numerical_input_sha256(problem) != original
    paths=getproperty.(implementation_record(Formulation()).sources, :path)
    @test all(path -> path in paths,
        (
            "src/datamodel/baseparams/geometry.jl",
            "src/datamodel/baseparams/resistance.jl",
            "src/datamodel/baseparams/inductance.jl",
            "src/datamodel/baseparams/dielectrics.jl",
            "src/datamodel/placement/paths.jl"))
    @test "src/engine/earthimpedance/formulas/carson1926.jl" in paths
    source=implementation_record()
    @test all(record -> bytes2hex(GauntletSupport.Gauntlet.sha256(record.source)) == record.sha256,source)
end
