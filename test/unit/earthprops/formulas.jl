@testitem "EarthProps / formula discovery and construction" tags=[:unit] begin
    const EP=LineCableModels.EarthProps
    const FD=EP.FD

    expected=(
        :AlipioVisacro2014,
        :CIGRE2019,
        :DatsiosMikropoulos2019,
        :LongmireSmith1975,
        :Messier1985,
        :Portela1999,
        :Scott1967,
        :VisacroAlipio2012,
        :VisacroPortela1987
    )
    @test FD.formulas() == expected
    files=sort(filter(endswith(".jl"), readdir(joinpath(
        pkgdir(LineCableModels), "src", "earthprops", "fd", "formulas"
    ))))
    @test files == collect(lowercase.(string.(expected)) .* ".jl")

    for identifier in expected
        formula=FD.Formula(identifier)
        @test formula_id(formula) === identifier
        @test FD.assumptions(formula) == FD.assumptions(Val(identifier))
        @test !isempty(description(formula))
        @test isconcretetype(typeof(formula))
    end
    @test_throws ArgumentError FD.Formula(:Unknown)
    @test_throws ArgumentError FD.Formula(:Portela1999; unknown = 1)
    @test FD.DEFAULT === nothing
    @test FD.Formula(:default) === nothing

    route=(material, frequency, values)->EP.EarthMaterial(
        material.rho / values.scale,
        material.eps_r + frequency / frequency,
        material.mu_r
    )
    experimental=FD.Formula(:Experiment, route, (scale = 2.0,))
    material=EP.EarthMaterial(100.0, 10.0, 1.0)
    @test constitutive(FD.Formula(:default), material, 50.0) === material
    @test_throws ArgumentError FD.Formula(:default; route)
    @test_throws ArgumentError FD.Formula(:default, route)
    output=@inferred constitutive(experimental, material, 50.0)
    @test (output.rho, output.eps_r, output.mu_r) == (50.0, 11.0, 1.0)
end

@testitem "EarthProps / legacy soil laws / numerical references" tags=[:unit] begin
    using TOML
    const EP=LineCableModels.EarthProps
    const FD=EP.FD

    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "earth_properties.toml"
    ))
    material=EP.EarthMaterial(
        reference["input_rho_ohm_m"],
        reference["input_eps_r"],
        reference["input_mu_r"]
    )
    sample_frequencies=reference["frequencies_hz"]

    for identifier in FD.formulas()
        formula=FD.Formula(identifier)
        expected=reference["formula"][String(identifier)]
        for (frequency, rho, eps_r) in zip(
            sample_frequencies,
            expected["rho_ohm_m"],
            expected["eps_r"]
        )
            output=@inferred constitutive(formula, material, frequency)
            @test output.rho ≈ rho rtol=2e-9
            @test output.eps_r ≈ eps_r rtol=2e-9
            @test output.mu_r === material.mu_r
        end
    end

    bounded=constitutive(FD.Formula(:VisacroAlipio2012), material, 50.0)
    at_reference=constitutive(FD.Formula(:VisacroAlipio2012), material, 100.0)
    @test bounded == at_reference
end

@testitem "EarthProps / constitutive precision and uncertainty" tags=[:unit] begin
    using Measurements: Measurement, measurement, uncertainty
    const EP=LineCableModels.EarthProps
    const FD=EP.FD

    for identifier in FD.formulas()
        formula=FD.Formula(identifier)
        material32=EP.EarthMaterial(100.0f0, 10.0f0, 1.0f0)
        output32=@inferred constitutive(formula, material32, 1000.0f0)
        @test output32 isa EP.EarthMaterial{Float32}

        uncertain_material=EP.EarthMaterial(
            measurement(100.0, 2.0),
            measurement(10.0, 0.5),
            measurement(1.0, 0.01)
        )
        uncertain_frequency=measurement(1000.0, 1.0)
        output=@inferred constitutive(formula, uncertain_material, uncertain_frequency)
        @test output isa EP.EarthMaterial{Measurement{Float64}}
        @test uncertainty(output.rho) > 0
        @test uncertainty(output.eps_r) > 0
        @test uncertainty(output.mu_r) == 0.01
    end
end

@testitem "Engine / earth constitutive cutover and air pass-through" tags=[:unit] setup=[
    UseEngineSupport,
] begin
    model=EarthModel(100.0, 10.0, 1.0)
    frequencies=[100.0, 1.0e6]
    formulation=Formulation(earth_properties = :CIGRE2019)
    @test formula_id(formulation.methods.earth_properties) === :CIGRE2019
    evaluated=@inferred LineCableModels.Engine._earth_data(
        formulation,
        (earth = model, freq = frequencies)
    )

    @test evaluated.rho[1, :] == fill(Inf, 2)
    @test evaluated.eps_r[1, :] == fill(1.0, 2)
    @test evaluated.mu_r[1, :] == fill(1.0, 2)
    @test evaluated.rho[2, 1] != evaluated.rho[2, 2]
    @test evaluated.eps_r[2, 1] != evaluated.eps_r[2, 2]
    @test evaluated.mu_r[2, :] == fill(1.0, 2)

    scan=collect(range(100.0, 1.0e6; length = 1000))
    scan_input=(earth = model, freq = scan)
    static_formulation=Formulation()
    LineCableModels.Engine._earth_data(static_formulation, scan_input)
    LineCableModels.Engine._earth_data(formulation, scan_input)
    static_allocation=@allocated LineCableModels.Engine._earth_data(
        static_formulation,
        scan_input
    )
    dispersive_allocation=@allocated LineCableModels.Engine._earth_data(
        formulation,
        scan_input
    )
    @test dispersive_allocation <= static_allocation + 1024
end
