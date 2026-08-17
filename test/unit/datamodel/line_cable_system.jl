@testitem "DataModel / LineCableSystem / derived counts and strict insertion" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    mapping(phase) = Dict("core"=>phase, "sheath"=>0, "jacket"=>0)
    first_position=CablePosition(design, 0.0, -1.0, mapping(1))
    system=LineCableSystem("line", 1000.0, first_position)

    @test system isa LineCableSystem{Float64}
    @test ncables(system) == 1
    @test nphases(system) == 1
    @test !hasproperty(system, :num_cables)
    @test !hasproperty(system, :num_phases)
    @test validate(system) === system

    second_position=CablePosition(design, 0.1, -1.0, mapping(2))
    @test add!(system, second_position) === system
    @test ncables(system) == 2
    @test nphases(system) == 2

    snapshot=deepcopy(system)
    overlap=CablePosition(design, 0.1001, -1.0, mapping(3))
    @test_throws DomainError add!(system, overlap)
    @test ncables(system) == ncables(snapshot)

    position32=convert(CablePosition{Float32}, CablePosition(design, 0.2, -1.0, mapping(3)))
    @test_throws ArgumentError add!(system, position32)
    @test ncables(system) == ncables(snapshot)

    system32=convert(LineCableSystem{Float32}, system)
    @test system32 isa LineCableSystem{Float32}
    @test ncables(system32) == ncables(system)
    @test nphases(system32) == nphases(system)
    @test convert(LineCableSystem{Float64}, system) === system

    @test_throws DomainError LineCableSystem("bad", 0.0, first_position)
    @test_throws ArgumentError LineCableSystem("", 1.0, first_position)
end
