@testitem "Validation / Tubular / representative rule classes" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    const V=LineCableModels.Validation
    const TubularType=LineCableModels.DataModel.Tubular
    material=LineCableModels.Materials.Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)

    valid=V.validate!(TubularType, 0.01, 0.02, material; temperature = 20.0)
    @test valid.r_in == 0.01
    @test valid.r_ex == 0.02
    @test valid.material_props === material

    invalid_cases=(
        (()->V.validate!(TubularType, -0.01, 0.02, material), ArgumentError),
        (()->V.validate!(TubularType, 0.02, 0.01, material), ArgumentError),
        (()->V.validate!(TubularType, 0.01, Inf, material), DomainError),
        (()->V.validate!(TubularType, 0.01, 0.02, "copper"), ArgumentError),
        (
            ()->V.validate!(TubularType, 0.01, 0.02, material; temperature = NaN),
            DomainError
        )
    )
    for (operation, exception_type) in invalid_cases
        @test_throws exception_type operation()
    end

    rules=V._rules(TubularType)
    @test any(rule -> rule isa V.Less && rule.a == :r_in && rule.b == :r_ex, rules)
    @test any(rule -> rule isa V.IsA{LineCableModels.Materials.Material}, rules)
end
