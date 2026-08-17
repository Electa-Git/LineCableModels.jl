@testitem "DataModel / CableDesign / strict insertion and reference state" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    source=TestFixtures.mv_cable_design()
    first_component=deepcopy(source.components[1])
    second_component=deepcopy(source.components[2])
    design=CableDesign("test", first_component)

    @test design isa CableDesign{Float64}
    @test design.nominal_data === nothing
    @test validate(design) === design
    @test add!(design, second_component) === design
    @test length(design.components) == 2
    @test design.components[2] === second_component

    snapshot=deepcopy(design)
    duplicate=CableComponent(
        second_component.id,
        deepcopy(second_component.conductor_group),
        deepcopy(second_component.insulator_group)
    )
    @test_throws ArgumentError add!(design, duplicate)
    @test length(design.components) == length(snapshot.components)

    component32=convert(CableComponent{Float32}, source.components[3])
    @test_throws ArgumentError add!(design, component32)
    @test length(design.components) == length(snapshot.components)

    bad_material=Material(1.0e-8, 1.0, 1.0, 25.0, 0.004)
    conductor=ConductorGroup(Tubular(second_component.insulator_group.r_ex, 0.04, bad_material))
    insulation=InsulatorGroup(Insulator(0.04, 0.045, Material(1e14, 2.3, 1.0, 25.0, 0.0)))
    @test_throws ArgumentError add!(design, CableComponent("hot", conductor, insulation))

    text=sprint(show, MIME("text/plain"), CableDesign("plain", first_component))
    @test occursin("CableDesign \"plain\"", text)
    @test !occursin("nominal values", text)

    design32=convert(CableDesign{Float32}, source)
    @test design32 isa CableDesign{Float32}
    @test convert(CableDesign{Float64}, source) === source
    @test length(design32.components) == length(source.components)
end
