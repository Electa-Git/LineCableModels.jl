@testitem "Earth / explicit reduction validates inputs and admits a user-owned rule" tags=[:unit] setup=[FormulaFixtures] begin
    const EP = LineCableModels.Earth
    const EH = EP.EquivalentHomogeneous
    const E = LineCableModels.Engine
    rule = EH.Formula(:default)
    @test EH.Formula(rule) === rule
    @test occursin("Bottommost", description(rule))
    @test occursin("before", description(EH.BeforeFD(:default)))
    @test occursin("after", description(EH.AfterFD(:default)))
    @test EH.rule(EH.BeforeFD(rule)) === rule
    @test EH.Formula(formula(:default)) isa typeof(rule)
    @test_throws ArgumentError EH.Formula(formula(:default; order = :before))
    @test_throws ArgumentError EH.Formula(formula(:default; equivalent_earth = formula(:default)))
    @test_throws ArgumentError EH.AbstractSequence(formula(:default; equivalent_earth = formula(:default)))
    @test_throws ArgumentError EH.Formula(:default;parameters=(rho=10.0,))
    model = build(EP.EarthModel, (EP.EarthLayer(100.0, 10.0, 1.0, 0.5),
        EP.EarthLayer(500.0, 20.0, 1.0)))
    pair = E.EarthPair(1, 2, (-0.25, -1.0), 0.75, (2, 3))
    rho = [Inf, 100.0, 500.0]
    epsilon = [1.0, 10.0, 20.0]
    mu = [1.0, 1.0, 1.0]
    @test_throws DimensionMismatch rule(rho[1:2], epsilon, mu, model, pair, 50.0)
    for frequency in (0.0, -50.0, Inf, NaN)
        @test_throws DomainError rule(rho, epsilon, mu, model, pair, frequency)
    end
    original = (copy(rho), copy(epsilon), copy(mu))
    selected=FormulaFixtures.MeanEarth()
    @test EH.Formula(selected) === selected
    material = selected(rho, epsilon, mu, model, pair, 50.0)
    @test length(selected.seen) == 1
    @test (material.rho, material.eps_r, material.mu_r) == (300.0, 15.0, 1.0)
    @test (rho, epsilon, mu) == original
    converted = convert(EP.EarthMaterial{Float32}, material)
    @test eltype(converted) === Float32
    @test eltype(typeof(converted)) === Float32
    @test (converted.rho, converted.eps_r, converted.mu_r) == (300.0f0, 15.0f0, 1.0f0)
    @test convert(EP.EarthMaterial{Float32}, converted) === converted
    restored = convert(EP.EarthMaterial{Float64}, converted)
    @test (restored.rho, restored.eps_r, restored.mu_r) == (material.rho, material.eps_r, material.mu_r)
end
