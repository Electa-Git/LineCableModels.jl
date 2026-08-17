@testitem "Validation / Tubular / owner-local declarative rules" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    const V=LineCableModels.Validation
    const TubularType=LineCableModels.DataModel.Tubular
    material=Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)

    layer=Tubular(0.01, 0.02, material)
    @test validate(layer) === layer
    @test V.check(TubularType, layer) === layer
    @test any(rule -> rule isa V.Less &&
                      rule.left == :r_in && rule.right == :r_ex,
        V.rules(TubularType))

    @test_throws DomainError Tubular(-0.01, 0.02, material)
    @test_throws DomainError Tubular(0.02, 0.01, material)
    @test_throws DomainError Tubular(0.01, Inf, material)
    @test_throws MethodError Tubular(0.01, 0.02, "copper")
    @test_throws MethodError Tubular(0.01, 0.02, material; temperature = NaN)
end
