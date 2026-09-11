@testitem "Materials / temperature law / equations, types and reference calibration" tags=[:unit] begin
    using Measurements
    const TD = LineCableModels.Materials.TemperatureDependent
    selected = TD.Formula(formula(:default))
    @test TD.formulas() == (:default,)
    @test formula_id(selected) === :default
    for T in (Float32, Float64, BigFloat)
        material = Material(:conductor, T(1.72e-8), one(T), one(T), T(20), T(0.004))
        for temperature in T.((20, 60, 80))
            rho = @inferred constitutive(selected, material, temperature)
            @test rho isa T
            @test rho ≈ material.rho * (1 + material.alpha * (temperature - material.T0))
            @test inv(rho) ≈ inv(material.rho) / (1 + material.alpha * (temperature - material.T0))
        end
        @test constitutive(nothing, material, T(250)) === material.rho
        @test material.rho == T(1.72e-8)
        @test material.T0 == T(20)
        @test material.alpha == T(0.004)
        @test_throws DomainError constitutive(selected, material, T(170))
        @test_throws DomainError constitutive(selected, material, T(-130))
        @test constitutive(selected, material, T(169)) > material.rho
        @test_throws DomainError constitutive(nothing, material, T(Inf))
    end
    fixed = Material(:conductor, 1.72e-8)
    @test constitutive(selected, fixed, 80.0) == fixed.rho
    passive = Material(:insulator, Inf, 2.3, 1, 20, -0.003)
    @test constitutive(selected, passive, 80.0) == Inf
    invalid = Material(:semicon, 1e3, 10, 1, 20, -0.05)
    @test_throws DomainError constitutive(selected, invalid, 40.0)
    uncertain_material = Material(:conductor, measurement(1.72e-8, 1e-10),
        1, 1, 20, 0.004)
    rho = constitutive(selected, uncertain_material, measurement(80.0, 1.0))
    @test nominal(rho) ≈ 1.72e-8 * 1.24
    @test uncertainty(rho) ≈ hypot(1.24e-10, 1.72e-8 * 0.004)
    @test uncertainty(uncertain_material.rho) == 1e-10
end

@testitem "Materials / temperature law / custom contributions and applicability" tags=[:unit] begin
    const TD = LineCableModels.Materials.TemperatureDependent
    const Binding = LineCableModels.FormulaMethod{:default, typeof(TD.temperature_resistivity)}
    calls = Ref(0)
    law = (material, temperature, parameters, options, workspace) -> begin
        calls[] += 1
        @test isempty(parameters) && isempty(options)
        @test workspace === nothing
        material.rho * exp((temperature - material.T0) / 1000)
    end
    @test_throws ArgumentError TD.Formula(:default; hooks=(contribution=law,))
    @eval LineCableModels.computation_options(::$Binding, ::$(typeof(law))) = (;)
    selected = TD.Formula(formula(:default; hooks=(contribution=law,)))
    material = Material(:conductor, 1.72e-8, 1, 1, 20, 0.004)
    @test constitutive(selected, material, 250.0) ≈ material.rho * exp(0.23)
    @test calls[] == 1
    @test selected.hooks.contribution === law
    @test_throws ArgumentError TD.Formula(:unknown)
    @test_throws ArgumentError TD.Formula(:default; parameters=(alpha=0.1,))
    @test_throws ArgumentError TD.Formula(:default; options=(integration=(method=:quad,),))
    @test_throws ArgumentError TD.Formula(:default; hooks=(other=law,))
    @test_throws ArgumentError TD.Formula(formula(:default; equivalent_earth=:default))
    for value in (0.0, -1.0, NaN, Inf, 1+im, [1.0])
        bad = (m, t, p, o, w) -> value
        @eval LineCableModels.computation_options(::$Binding, ::$(typeof(bad))) = (;)
        invalid = TD.Formula(:default; hooks=(contribution=bad,))
        @test_throws DomainError constitutive(invalid, material, 80.0)
    end
    for constructor in (Formulation, CableConstantsFormulation, LineCableModelsFEM)
        @test constructor(temperature_dependence=nothing).methods.temperature_dependence === nothing
        @test_throws ArgumentError constructor(options=(temperature_correction=true,))
        grid = constructor(temperature_dependence=Grid((formula(:default), nothing)))
        @test length(grid) == 2
        @test first(grid).methods.temperature_dependence isa TD.Formula
        @test last(collect(grid)).methods.temperature_dependence === nothing
    end
end
