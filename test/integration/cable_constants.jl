@testitem "DataModel / cable constants / datasheet and physical invariants" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    design=TestFixtures.mv_cable_design()
    constants=compute(CableConstantsProblem(design), Formulation())
    nominal=design.nominal_data

    @test constants.R > 0
    @test constants.L > 0
    @test constants.C > 0
    @test basis(constants) === :per_length
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C

    @test isapprox(1_000 * constants.R, nominal.resistance; rtol = 0.06)
    @test isapprox(1e6 * constants.L, nominal.inductance; rtol = 0.06)
    @test isapprox(1e9 * constants.C, nominal.capacitance; rtol = 0.06)

    separated=compute(
        CableConstantsProblem(design; separation = 2.5*outer_radius(design)),
        Formulation()
    )
    @test separated.R == constants.R
    @test separated.C == constants.C
    @test separated.L > constants.L
end

@testitem "Fixtures / factories / mutable state is never shared" tags=[:integration] setup=[
    TestFixtures,
] begin
    first_design=TestFixtures.mv_cable_design()
    second_design=TestFixtures.mv_cable_design()
    @test first_design !== second_design
    @test first_design.components !== second_design.components
    @test first_design.components[1] !== second_design.components[1]

    first_system=TestFixtures.three_phase_system()
    second_system=TestFixtures.three_phase_system()
    @test first_system !== second_system
    @test first_system.cables !== second_system.cables
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
