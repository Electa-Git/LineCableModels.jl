@testitem "DataModel / InsulatorGroup / reference-frequency insertion" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    first_layer=Insulator(0.006, 0.008, insulator_props)
    second_layer=Insulator(0.008, 0.010, insulator_props)
    group=InsulatorGroup(first_layer; reference_frequency = 60.0)
    ω=2π*group.reference_frequency
    expected=parallel(
        complex(group.shunt_conductance, ω*group.shunt_capacitance),
        complex(second_layer.shunt_conductance, ω*second_layer.shunt_capacitance)
    )

    @test add!(group, second_layer) === group
    @test length(group.layers) == 2
    @test group.r_ex == second_layer.r_ex
    @test group.shunt_conductance ≈ real(expected)
    @test group.shunt_capacitance ≈ imag(expected)/ω
    @test group.reference_frequency == 60.0
    @test validate(group) === group

    snapshot=deepcopy(group)
    @test_throws DomainError add!(group, Insulator(0.011, 0.012, insulator_props))
    @test group.r_ex == snapshot.r_ex
    @test length(group.layers) == length(snapshot.layers)

    material32=convert(Material{Float32}, insulator_props)
    layer32=Insulator(Float32(0.010), Float32(0.011), material32)
    @test_throws ArgumentError add!(group, layer32)
    @test group.r_ex == snapshot.r_ex

    converted=convert(InsulatorGroup{Float32}, group)
    @test converted isa InsulatorGroup{Float32}
    @test converted.reference_frequency isa Float32
    @test convert(InsulatorGroup{Float64}, group) === group
end
