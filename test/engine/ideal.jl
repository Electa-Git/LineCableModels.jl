# The images of an earth of infinite conductivity on the overhead reference. Its two wires of
# radius 42.5 mm stand 1 m above the interface and 1 m apart.
@testset "Engine ideal earth" begin
    Engine = LineCableModels.Engine
    overhead = Dict(REFERENCES)["overhead_bare_wires"]
    r, h, D = 42.5e-3, 1.0, 1.0
    ε₀ = 8.8541878128e-12
    ideal = compute(overhead, CoaxialFormulation(earth_admittance = :ideal); trace = true)
    both = compute(overhead, CoaxialFormulation(); trace = true)
    @test @inferred(compute(overhead, CoaxialFormulation(earth_admittance = :ideal))) isa
          Engine.LineParameters
    @test isequal(ideal.Z, both.Z)
    Pg = ideal.trace.Pg
    @test all(isapprox(Pg[i, i, k], log(2h / r) / (2π * ε₀); rtol = 1e-15)
    for i in 1:2, k in axes(Pg, 3))
    @test all(isapprox(Pg[1, 2, k], log(hypot(D, 2h) / D) / (2π * ε₀); rtol = 1e-15)
    for k in axes(Pg, 3))
    @test_throws "gives no impedance!" CoaxialFormulation(earth_impedance = :ideal)
    # The earth shields the buried wire of the mixed pair.
    mixed = Dict(REFERENCES)["mixed_bare_wires"]
    input = Engine.CoaxialInput(mixed, CoaxialFormulation().options)
    formula = Engine.EarthReturn.Formula(Val(:ideal))
    s = im * (2π * 50.0)
    @test iszero(formula(Val(:P), Val(:self), 2, (;), input, mixed.earth.layers, s, 1))
    @test iszero(formula(Val(:P), Val(:mutual), 1, 2, (;), input, mixed.earth.layers, s, 1))
end
