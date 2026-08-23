@testitem "ParametricBuilder / SystemBuilder / flat formations and phase maps" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    const PB=LineCableModels.ParametricBuilder

    design=TestFixtures.mv_cable_design()
    horizontal=PB.hflat(
        x = 0.0,
        y = -1.0,
        spacing = 0.10,
        phases = [:core=>(1, 2)],
        n = 2
    )
    vertical=PB.vflat(
        x = 1.0,
        y = -1.0,
        spacing = 0.10,
        phases = (:core=>[3, 4],),
        n = 2
    )
    earth=EarthModel(100.0, 10.0, 1.0)
    specification=PB.SystemBuilder(
        "flat-formations",
        design,
        (horizontal, [vertical]);
        earth,
        frequencies = [50.0]
    )
    problem=only(specification)
    positions=problem.system.cables

    @test length(positions) == 4
    @test [(position.horz, position.vert) for position in positions] ==
          [(0.0, -1.0), (0.1, -1.0), (1.0, -1.0), (1.0, -1.1)]
    @test map(position -> first(position.conn), positions) == [1, 2, 3, 4]
    @test problem.earth_props.layers[end].rho == earth.layers[end].rho

    unconnected=only(PB.at(x = 0.0, y = -1.0, phases = nothing))
    @test only(unconnected.connections) == Dict{String, Int}()
    uniform=only(PB.at(x = 0.0, y = -1.0, phases = :core=>2))
    @test only(uniform.connections) == Dict("core" => 2)

    @test_throws ArgumentError PB.hflat(
        spacing = 0.1,
        phases = :core => (),
        n = 0
    )
    @test_throws ArgumentError PB.vflat(
        spacing = 0.1,
        phases = :core => (),
        n = 0
    )
    @test_throws DimensionMismatch PB.hflat(
        spacing = 0.1,
        phases = :core => (1,),
        n = 2
    )
    @test_throws ArgumentError PB.at(x = 0.0, y = -1.0, phases = :core => 1.5)
    @test_throws ArgumentError PB.PositionDefinition(
        Val(:point),
        ("not-real", -1.0),
        (Dict("core" => 1),)
    )

    too_close=only(PB.hflat(
        spacing = eps(),
        phases = :core=>(1, 2),
        n = 2
    ))
    @test_throws ArgumentError PB._position_coordinates(too_close, design)
    unsupported=PB.PositionDefinition(
        Val(:unsupported),
        (0.0, -1.0, 0.10),
        (Dict("core"=>1),)
    )
    @test_throws ArgumentError PB._position_coordinates(unsupported, design)
    @test_throws ArgumentError PB.SystemBuilder(
        "invalid-positions",
        design,
        42;
        earth,
        frequencies = [50.0]
    )
end
