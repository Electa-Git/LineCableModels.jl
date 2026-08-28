@testitem "DataModel / LineCableSystem / completed immutable state" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    mapping(phase) = Dict(
        terminal=>(index==1 ? phase : 0)
    for (index, terminal) in enumerate(design.terminal_order)
    )
    system=build(
        LineCableSystem,
        [design, design],
        [Pose2(0.0, -1.0, 0.0), Pose2(0.1, -1.0, 0.0)];
        connections = [mapping(1), mapping(2)],
        system_id = "line",
        line_length = 1000.0
    )

    @test system isa LineCableSystem{Float64}
    @test ncables(system) == 2
    @test nphases(system) == 2
    @test system.designs == [design, design]
    @test system.positions == [Pose2(0.0, -1.0, 0.0), Pose2(0.1, -1.0, 0.0)]
    @test system.connection_order == vcat(
        [mapping(1)[terminal] for terminal in design.terminal_order],
        [mapping(2)[terminal] for terminal in design.terminal_order]
    )
    @test length(system.geometry) == 2 * length(design.geometry.regions)
    @test !hasproperty(system, :num_cables)
    @test !hasproperty(system, :num_phases)
    @test validate(system) === system

    @test_throws DomainError build(
        LineCableSystem,
        [design, design],
        [Pose2(0.0, -1.0, 0.0), Pose2(0.001, -1.0, 0.0)];
        connections = [mapping(1), mapping(2)]
    )
    @test_throws DimensionMismatch build(
        LineCableSystem,
        [design, design],
        [Pose2(0.0, -1.0, 0.0)]
    )

    @test_throws DomainError build(
        LineCableSystem, design, Pose2(0.0, -1.0, 0.0); line_length = 0.0
    )
    @test_throws ArgumentError build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0, 0.0);
        system_id = ""
    )
end
