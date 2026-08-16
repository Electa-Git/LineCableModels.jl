@testitem "DataModel / presentation / tree summaries and tabular mappings" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    system=TestFixtures.three_phase_system()
    component=first(design.components)
    conductor_part=first(component.conductor_group.layers)
    insulator_part=first(component.insulator_group.layers)

    @test eltype(design) === Float64
    @test eltype(typeof(design)) === Float64
    @test eltype(component) === Float64
    @test eltype(typeof(component)) === Float64
    @test eltype(first(system.cables)) === Float64
    @test eltype(system) === Float64

    component_text=sprint(show, MIME"text/plain"(), component)
    @test occursin("CableComponent \"$(component.id)\"", component_text)
    @test occursin("ConductorGroup", component_text)
    @test occursin("InsulatorGroup", component_text)
    @test occursin("Effective properties", component_text)

    design_text=sprint(show, MIME"text/plain"(), design)
    @test occursin("CableDesign \"$(design.cable_id)\"", design_text)
    @test all(part -> occursin("Component \"$(part.id)\"", design_text), design.components)
    @test occursin("nominal values", design_text)

    for object in (
        conductor_part,
        insulator_part,
        component.conductor_group,
        component.insulator_group
    )
        rendered=sprint(show, MIME"text/plain"(), object)
        @test occursin(string(nameof(typeof(object))), rendered)
        @test occursin("r_in=", rendered)
    end

    system_text=sprint(show, MIME"text/plain"(), system)
    @test occursin("LineCableSystem \"$(system.system_id)\"", system_text)
    @test occursin("core→1", system_text)
    @test occursin("core→3", system_text)

    system_table=DataFrame(system)
    @test nrow(system_table) == 3
    @test names(system_table) == ["cable_id", "horz", "vert", "phase_mapping"]
    @test system_table.phase_mapping == [
        "core: 1, sheath: 0, jacket: 0",
        "core: 2, sheath: 0, jacket: 0",
        "core: 3, sheath: 0, jacket: 0"
    ]

    library=@test_logs (:info, r"Initializing empty cables") CablesLibrary()
    library[design.cable_id]=design
    @test length(library) == 1
    @test first(values(library)) === design
    @test first(library)[1] == design.cable_id
    @test_logs (:info, r"loaded from") begin
        @test get(library, design.cable_id) === design
    end
    sentinel=Ref(:missing)
    @test_logs (:warn, r"not found") begin
        @test get(library, "missing", sentinel) === sentinel
    end

    library_table=DataFrame(library)
    @test library_table.cable_id == [design.cable_id]
    @test library_table.components == [join(getproperty.(design.components, :id), ", ")]
    @test_logs (:info, r"removed") delete!(library, design.cable_id)
    @test isempty(library)
    @test_logs (:error, r"cannot delete") @test_throws KeyError delete!(library, "missing")

    routed_part=RectStrands(
        0.01,
        0.001,
        0.002,
        12,
        12.0,
        TestFixtures.copper_material()
    )
    @test hasproperty(routed_part, :shape)
    @test hasproperty(routed_part, :width)
    @test :width in propertynames(routed_part)
    @test routed_part.width == routed_part.shape.width
    @test_throws FieldError routed_part.property_that_does_not_exist

    buffer=IOBuffer()
    displayed=LineCableModels.DataModel._print_fields(
        buffer,
        (; number = 1.23456, label = "fixture", ignored = NaN),
        [:number, :label, :ignored, :absent];
        sigdigits = 3
    )
    @test displayed == 2
    @test String(take!(buffer)) == "number=1.23, label=fixture"
end

@testitem "DataModel / simplification / geometry-preserving reconstruction" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    original=TestFixtures.mv_cable_design()
    equivalent_design=equivalent(original; new_id = "equivalent")
    simplified=LineCableModels.DataModel.nonsensify(original; new_id = "simplified")

    for rebuilt in (equivalent_design, simplified)
        @test rebuilt !== original
        @test length(rebuilt.components) == length(original.components)
        @test getproperty.(rebuilt.components, :id) ==
              getproperty.(original.components, :id)
        for (new_component, old_component) in zip(rebuilt.components, original.components)
            @test new_component.conductor_group.r_in == old_component.conductor_group.r_in
            @test new_component.conductor_group.r_ex == old_component.conductor_group.r_ex
            @test new_component.insulator_group.r_ex == old_component.insulator_group.r_ex
            @test length(new_component.conductor_group.layers) == 1
            @test length(new_component.insulator_group.layers) == 1
        end
    end

    @test equivalent_design.cable_id == "equivalent"
    @test simplified.cable_id == "simplified"
    @test equivalent(original).cable_id == "$(original.cable_id)_equivalent"
    @test LineCableModels.DataModel.nonsensify(original).cable_id ==
          "$(original.cable_id)_nonsense"

    empty_design=CableDesign{Float64}(
        "empty",
        CableComponent{Float64}[];
        nominal_data = NominalData()
    )
    @test_throws ArgumentError equivalent(empty_design)
    @test_throws ArgumentError LineCableModels.DataModel.nonsensify(empty_design)
end
