@testitem "DataModel / LineCableSystem / completed immutable state" tags=[:unit, :parametric] setup=[
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.coaxial_design()
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
    @test validate(system) === system

    neutral=build(
        LineCableSystem,
        design,
        Pose2(0.0, 0.0, 0.0);
        environment = nothing,
        system_id = "neutral"
    )
    @test only(neutral.positions) == Pose2(0.0, 0.0, 0.0)
    @test neutral.environment === nothing
    @test_throws OverlapError build(
        LineCableSystem,
        design,
        Pose2(0.0, 0.0, 0.0);
        environment = :interface
    )

    @test_throws OverlapError build(
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

@testitem "DataModel / LineCableSystem / untyped collections build by element" tags=[:unit, :parametric] setup=[
    UseDataModelSupport,
    TestFixtures
] begin
    design = TestFixtures.coaxial_design()
    poses = [Pose2(0.0, -1.0, 0.0), Pose2(0.1, -1.0, 0.0)]
    mapping(phase, T) = Dict{Symbol, T}(
        terminal => (index == 1 ? phase : 0)
    for (index, terminal) in enumerate(design.terminal_order)
    )
    typed = build(LineCableSystem, [design, design], poses;
        connections = [mapping(1, Int), mapping(2, Int)])
    # A Vector{Any} of designs and a Dict{Symbol, Any} of integers build the same system.
    designs = []
    push!(designs, design, design)
    untyped = build(LineCableSystem, designs, Any[poses...];
        connections = Any[mapping(1, Any), mapping(2, Any)])
    @test untyped.designs == typed.designs
    @test untyped.positions == typed.positions
    @test untyped.connections == typed.connections
    @test untyped.connection_order == typed.connection_order
    @test untyped.geometry == typed.geometry
    # An offending element is named by its type.
    @test_throws "designs must contain completed CableDesign objects; got Int64" build(
        LineCableSystem, Any[design, 1], poses)
    message = "placements must contain Pose2 values or `(x, y[, φ])` tuples; got String"
    @test_throws message build(LineCableSystem, [design, design], Any[first(poses), "east"])
    fractional = mapping(1, Any)
    fractional[first(design.terminal_order)] = 1.5
    @test_throws "connection assignments must be integers; got Float64" build(
        LineCableSystem, design, first(poses); connections = fractional)
end

@testitem "DataModel / LineCableSystem / terminal labels" tags=[:unit, :parametric] setup=[
    UseDataModelSupport,
    TestFixtures
] begin
    design = TestFixtures.coaxial_design()
    system = build(LineCableSystem, [design, design],
        [Pose2(0.0, -1.0, 0.0), Pose2(0.1, -1.0, 0.0)];
        connections = [Dict(:core => 7, :sheath => 0), Dict(:core => 7, :sheath => 3)])

    @test terminal_labels(system) ==
        ["cable:1:core", "cable:1:sheath", "cable:2:core", "cable:2:sheath"]
    @test terminal_labels(system, [4, 1, 2]; bundle = false) ==
        ["cable:2:sheath", "cable:1:core", "cable:1:sheath"]
    @test terminal_labels(system, [4, 1, 2]; bundle = true) ==
        ["cable:2:sheath", "bundle:[cable:1:core,cable:2:core]", "cable:1:sheath"]
    @test terminal_labels(system, [3]; bundle = true) ==
        ["bundle:[cable:1:core,cable:2:core]"]
    separate = build(LineCableSystem, [design, design], system.positions;
        connections = [Dict(:core => 7, :sheath => 0), Dict(:core => 3, :sheath => 0)])
    @test terminal_labels(separate, [2, 3, 4, 1]; bundle = true) ==
        ["cable:1:sheath", "cable:2:core", "cable:2:sheath", "cable:1:core"]
    @test terminal_labels(system, Int[]; bundle = true) == String[]
    @test terminal_labels(system, Int[]; bundle = false) == String[]
    @test terminal_labels((cable = cable, terminal = :core) for cable in (2, 1)) ==
        ["cable:2:core", "cable:1:core"]
end
