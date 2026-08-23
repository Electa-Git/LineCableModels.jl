@testitem "ParametricBuilder / CableBuilder / wire, strip, and semiconductor staging" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    import LineCableModels.ParametricBuilder as PB

    conductor_material=TestFixtures.copper_material()
    dielectric_material=TestFixtures.insulator_material()
    semiconductor_material=TestFixtures.semicon_material()

    stranded=PB.Conductor.Stranded(
        :core;
        layers = 3,
        wire_radius = 0.001,
        num_wires = 6,
        lay_ratio = 11.0,
        material = conductor_material
    )
    @test length(stranded) == 2
    @test length(PB.Conductor.Stranded(
        :single;
        layers = 1,
        wire_radius = 0.001,
        material = conductor_material
    )) == 1
    core_insulation=PB.Insulator.Tubular(
        :core;
        thickness = 0.004,
        material = dielectric_material
    )
    screen=PB.Conductor.Strip(
        :screen;
        thickness = 0.001,
        width = 0.0005,
        lay_ratio = 10.0,
        material = conductor_material
    )
    screen_semicon=PB.Insulator.Semicon(
        :screen;
        thickness = 0.001,
        material = semiconductor_material
    )
    design=only(PB.Gridspace(PB.CableBuilder(
        "staged-wire-strip",
        stranded,
        core_insulation,
        [screen, screen_semicon],
        nominal = NominalData()
    )))

    @test design.cable_id == "staged-wire-strip"
    @test length(design.components) == 2
    @test length(design.components[1].conductor_group.layers) == 3
    @test design.components[1].conductor_group.layers[1] isa CircStrands
    @test design.components[1].conductor_group.layers[1].num_wires == 1
    @test design.components[1].conductor_group.layers[2].num_wires == 6
    @test design.components[1].conductor_group.layers[3].num_wires == 12
    @test only(design.components[2].conductor_group.layers) isa
          LineCableModels.DataModel.Strip
    @test only(design.components[2].insulator_group.layers) isa Semicon
    @test design.components[2].conductor_group.r_in ==
          design.components[1].insulator_group.r_ex

    named_nominal=PB._nominal_data((designation_code = "typed", U0 = 12.0))
    @test named_nominal isa NominalData
    @test named_nominal.designation_code == "typed"
    @test PB._nominal_data(named_nominal) === named_nominal
    @test PB._nominal_data(nothing) === nothing
end

@testitem "ParametricBuilder / CableBuilder / validation equivalence classes" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    import LineCableModels.ParametricBuilder as PB

    conductor_material=TestFixtures.copper_material()
    dielectric_material=TestFixtures.insulator_material()

    wire_cases=(
        (num_wires = 0, lay_ratio = 10.0),
        (num_wires = 6, lay_ratio = -1.0)
    )
    for case in wire_cases
        grid=PB.Conductor.Wires(
            :core;
            wire_radius = 0.001,
            num_wires = case.num_wires,
            lay_ratio = case.lay_ratio,
            material = conductor_material
        )
        @test_throws ArgumentError first(grid)
    end

    strip_cases=(
        (width = 0.0, lay_ratio = 10.0),
        (width = 0.001, lay_ratio = -1.0)
    )
    for case in strip_cases
        grid=PB.Conductor.Strip(
            :screen;
            thickness = 0.001,
            width = case.width,
            lay_ratio = case.lay_ratio,
            material = conductor_material
        )
        @test_throws ArgumentError first(grid)
    end

    @test_throws ArgumentError PB.Conductor.Stranded(
        :core;
        layers = 0,
        wire_radius = 0.001,
        material = conductor_material
    )
    @test_throws ArgumentError PB.Conductor.Tubular(
        :core;
        radius = 0.01,
        thickness = 0.001,
        material = conductor_material
    )
    @test_throws ArgumentError PB.Conductor.Tubular(
        :core;
        material = conductor_material
    )
    @test_throws ArgumentError PB.CableBuilder("empty")
    @test_throws ArgumentError PB.CableBuilder("wrong-part", 3)
    @test_throws ArgumentError PB.CableBuilder("bad-combine"; combine = :invalid)
    @test_throws ArgumentError PB._nominal_data(3)

    conductor=PB.Conductor.Solid(
        :core;
        radius = 0.01,
        material = conductor_material
    )
    insulation=PB.Insulator.Tubular(
        :core;
        thickness = 0.002,
        material = dielectric_material
    )
    @test_throws ArgumentError only(PB.Gridspace(PB.CableBuilder(
        "bad-nominal",
        conductor,
        insulation;
        nominal = 3
    )))
    @test_throws ArgumentError only(PB.Gridspace(PB.CableBuilder("no-insulator", conductor)))
    @test_throws ArgumentError only(PB.Gridspace(PB.CableBuilder("no-conductor", insulation)))

    @test_throws ArgumentError PB.PartDefinition(
        Val(:conductor),
        Val(LineCableModels.DataModel.Tubular),
        Val(:unsupported),
        :core,
        1,
        0.01,
        (),
        conductor_material
    )

    bad_part=PB.PartDefinition(
        Val(:conductor),
        Val(Int),
        Val(:radius),
        :core,
        1,
        0.01,
        (),
        conductor_material
    )
    @test_throws ArgumentError PB._materialize_part(bad_part, 0.0, 1)
end
