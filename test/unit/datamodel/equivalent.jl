@testitem "DataModel / equivalent and nonsensify / geometry-preserving reconstruction" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    source=TestFixtures.mv_cable_design()
    flattened=equivalent(source)
    simplified=LineCableModels.DataModel.nonsensify(source)

    @test flattened.cable_id == source.cable_id*"_equivalent"
    @test simplified.cable_id == source.cable_id*"_nonsense"
    @test length(flattened.components) == length(source.components)
    @test length(simplified.components) == length(source.components)
    @test flattened.nominal_data === source.nominal_data
    @test simplified.nominal_data === source.nominal_data
    @test validate(flattened) === flattened
    @test validate(simplified) === simplified

    for (original, reduced, plain) in zip(
        source.components,
        flattened.components,
        simplified.components
    )
        @test length(reduced.conductor_group.layers) == 1
        @test length(reduced.insulator_group.layers) == 1
        @test length(plain.conductor_group.layers) == 1
        @test length(plain.insulator_group.layers) == 1
        @test reduced.conductor_group.r_in == original.conductor_group.r_in
        @test reduced.conductor_group.r_ex == original.conductor_group.r_ex
        @test reduced.insulator_group.r_ex == original.insulator_group.r_ex
        @test reduced.conductor_props.rho ≈ original.conductor_props.rho
        @test reduced.insulator_props.eps_r ≈ original.insulator_props.eps_r
        @test plain.conductor_group.r_in == original.conductor_group.r_in
        @test plain.insulator_group.r_ex == original.insulator_group.r_ex
    end

    @test equivalent(source; new_id = "flat").cable_id == "flat"
    @test LineCableModels.DataModel.nonsensify(source; new_id = "plain").cable_id == "plain"
end
