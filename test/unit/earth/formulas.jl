@testitem "Earth / static relation and explicit constitutive hooks" tags=[:unit] begin
    using Measurements: measurement, uncertainty
    const EP=LineCableModels.Earth
    relation=EP.FrequencyDependent.Formula(:default)
    for T in (Float32, Float64, BigFloat)
        material=EP.EarthMaterial(T(100), T(10), T(1))
        @test relation(material, T(50)) === material
    end
    material=EP.EarthMaterial(measurement(100.0, 2.0), measurement(10.0, 0.2), measurement(1.0, 0.01))
    @test uncertainty(relation(material, 50.0).rho) == 2.0
    seen=Float64[]
    law=(material,
        frequency,
        parameters, options,
        workspace)->begin
        push!(seen, frequency)
        EP.EarthMaterial(material.rho/(1+frequency/100), material.eps_r, material.mu_r)
    end
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{
            :default, typeof(EP.FrequencyDependent.earth_material)},
        ::$(typeof(law))) = (;)
    selected=EP.FrequencyDependent.Formula(formula(:default; hooks = (contribution = law,)))
    @test selected.hooks.contribution === law
    @test selected(EP.EarthMaterial(100.0, 10.0, 1.0), 100.0).rho == 50.0
    @test seen == [100.0]
    @test_throws DomainError selected(EP.EarthMaterial(100.0, 10.0, 1.0), 0.0)
    @test_throws ArgumentError EP.FrequencyDependent.Formula(:default; parameters = (unknown = 1,))
end
