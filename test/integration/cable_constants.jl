@testitem "DataModel / cable constants / datasheet and physical invariants" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    design=TestFixtures.mv_cable_design()
    constants=CableConstants(design)

    @test constants.R > 0
    @test constants.L > 0
    @test constants.C > 0
    @test basis(constants) === :pul
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C

    altered=CableConstants(
        design;
        position = at(x = 0, y = -2),
        earth_props = Earth(rho = 200)
    )
    @test altered.R != constants.R
    @test altered.C ≈ constants.C
    @test altered.L != constants.L

    designs=build(
        CableDesign,
        Grid(("constants-a", "constants-b")),
        design.root
    )
    constant_space=CableConstants(designs)
    @test constant_space isa Gridspace{CableConstants}
    @test length(constant_space) == 2
    @test all(value -> value isa CableConstants, constant_space)
    @test first(constant_space) == CableConstants(first(designs))
end

@testitem "Fixtures / factories / mutable state is never shared" tags=[:integration] setup=[
    TestFixtures,
] begin
    first_design=TestFixtures.mv_cable_design()
    second_design=TestFixtures.mv_cable_design()
    @test first_design !== second_design
    @test first_design.root !== second_design.root
    @test first_design.geometry.regions !== second_design.geometry.regions

    first_system=TestFixtures.three_phase_system()
    second_system=TestFixtures.three_phase_system()
    @test first_system !== second_system
    @test first_system.designs !== second_system.designs
    @test first_system.positions !== second_system.positions
    @test ncables(first_system) == ncables(second_system) == 3
    @test nphases(first_system) == nphases(second_system) == 3

    first_parameters=TestFixtures.two_conductor_results()
    second_parameters=TestFixtures.two_conductor_results()
    @test first_parameters !== second_parameters
    @test first_parameters.Z.values !== second_parameters.Z.values
    first_parameters.Z.values[1]=complex(Inf, Inf)
    @test all(isfinite, second_parameters.Z.values)

    first_monte_carlo=TestFixtures.cable_monte_carlo_result()
    second_monte_carlo=TestFixtures.cable_monte_carlo_result()
    @test first_monte_carlo !== second_monte_carlo
    first_samples=only(samples(first_monte_carlo))
    second_samples=only(samples(second_monte_carlo))
    @test first_samples.R !== second_samples.R
    first_samples.R[1]=Inf
    @test all(isfinite, second_samples.R)
end
