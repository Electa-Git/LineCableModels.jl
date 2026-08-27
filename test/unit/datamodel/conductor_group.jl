@testitem "DataModel / ConductorGroup / strict atomic insertion" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    first_layer=Tubular(0.0, 0.005, copper_props)
    second_layer=Tubular(0.005, 0.006, copper_props)
    group=ConductorGroup(first_layer)
    expected_resistance=parallel(first_layer.resistance, second_layer.resistance)
    expected_alpha=equivalent_alpha(
        group.alpha,
        group.resistance,
        second_layer.material_props.alpha,
        second_layer.resistance
    )

    @test add!(group, second_layer) === group
    @test length(group.layers) == 2
    @test group.r_in == 0.0
    @test group.r_ex == 0.006
    @test group.resistance ≈ expected_resistance
    @test group.alpha ≈ expected_alpha
    @test validate(group) === group

    snapshot=deepcopy(group)
    @test_throws DomainError add!(group, Tubular(0.007, 0.008, copper_props))
    @test group.r_ex == snapshot.r_ex
    @test length(group.layers) == length(snapshot.layers)

    material32=convert(Material{Float32}, copper_props)
    layer32=Tubular(Float32(0.006), Float32(0.007), material32)
    @test_throws ArgumentError add!(group, layer32)
    @test group.r_ex == snapshot.r_ex

    different_reference=Material(:conductor, 1.7241e-8, 1.0, 1.0, 25.0, 0.00393)
    @test_throws ArgumentError add!(
        group,
        Tubular(group.r_ex, 0.007, different_reference)
    )

    converted=convert(ConductorGroup{Float32}, group)
    @test converted isa ConductorGroup{Float32}
    @test converted !== group
    @test convert(ConductorGroup{Float64}, group) === group
end
