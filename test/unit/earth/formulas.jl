@testitem "Earth / static and user-owned frequency-dependent relations" tags=[:unit, :earth] setup=[FormulaFixtures] begin
    using Measurements: measurement,uncertainty
    const EP=LineCableModels.Earth
    relation=EP.FrequencyDependent.Formula(:default)
    for T in (Float32,Float64,BigFloat)
        material=EP.EarthMaterial(T(100),T(10),T(1))
        @test relation(material,T(50)) === material
    end
    material=EP.EarthMaterial(measurement(100.0,2.0),measurement(10.0,0.2),measurement(1.0,0.01))
    @test uncertainty(relation(material,50.0).rho)==2.0
    selected=FormulaFixtures.DispersiveEarth()
    @test EP.FrequencyDependent.Formula(selected) === selected
    @test selected(EP.EarthMaterial(100.0,10.0,1.0),100.0).rho==50.0
    @test only(selected.seen)[2]==100.0
    @test_throws DomainError selected(EP.EarthMaterial(100.0,10.0,1.0),0.0)
    @test_throws ArgumentError EP.FrequencyDependent.Formula(:default;parameters=(unknown=1,))
end

@testitem "Earth / soil coefficients are checked at construction" tags=[:unit, :earth] begin
    const FD=LineCableModels.Earth.FrequencyDependent
    for id in FD.formulas(FD.Formula)
        selected=FD.Formula(id)
        for name in keys(selected.parameters), value in ("bad",true,NaN,Inf,1+im)
            @test_throws ArgumentError FD.Formula(id;parameters=NamedTuple{(name,)}((value,)))
        end
    end
    for (id,name,value) in ((:longmire1975,:corner_scale,0.0),
            (:visacro1987,:normalization_frequency,-1.0),
            (:visacro2012,:frequency_boundary,0.0),
            (:messier1985,:epsilon_infinity,-1.0),
            (:portela1999,:exponent,1.0),(:alipio2014,:exponent,-1.0))
        @test_throws ArgumentError FD.Formula(id;parameters=NamedTuple{(name,)}((value,)))
    end
end

@testitem "Earth / published dispersive relations match numerical references" tags=[:unit, :earth] begin
    const EP=LineCableModels.Earth
    const FD=EP.FrequencyDependent
    material=EP.EarthMaterial(100.0,10.0,1.0)
    # Independently evaluated from the published expressions documented by
    # each formula, for rho_0 = 100 ohm m and f = 10 kHz.
    references=(
        alipio2014=(98.085845947624506,409.88833451812962),
        cigre2019=(98.078364240070741,408.02591429681854),
        datsios2019=(99.953087658615345,77.228148054469685),
        longmire1975=(91.26463633360153,548.99524080259334),
        messier1985=(97.1029453119354,544.28505232978239),
        portela1999=(97.229937848620864,1088.2832591726246),
        scott1967=(101.39113857366789,542.00089040162334),
        visacro1987=(71.779429127136169,815.10931718426775),
        visacro2012=(98.649670481737814,192.20336879472805),
    )
    for (identifier,(rho,eps_r)) in pairs(references)
        evaluated=FD.Formula(identifier)(material,1.0e4)
        @test evaluated.rho ≈ rho rtol=2e-14
        @test evaluated.eps_r ≈ eps_r rtol=2e-14
        @test evaluated.mu_r === material.mu_r
    end

    visacro=FD.Formula(:visacro2012)
    @test visacro(material,50.0) == visacro(material,100.0)
    datsios=FD.Formula(:datsios2019)
    @test datsios(material,50.0).eps_r == datsios(material,3000.0).eps_r
end

@testitem "Earth / causal models and fitted real-frequency soil laws" tags=[:unit, :earth] begin
    const EP = LineCableModels.Earth
    const FD = EP.FrequencyDependent
    ε₀ = LineCableModels.Commons.vacuum_permittivity(Float64)
    coefficients = (3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1,
        1.25e1, 4.8, 2.17, 0.98, 0.392, 0.173)
    # Independent published real-frequency expressions.
    function published(identifier, rho, f)
        σ₀ = inv(rho)
        if identifier === :alipio2014
            σ̂₀ = 1000 / rho
            D = 1.26 * σ̂₀^(-0.73)
            σ = σ̂₀ * (1 + D * (f / 1e6)^0.54) / 1000
            εr = 12 + tan(π * 0.54 / 2) * 1e-3 * σ̂₀ * D * f^(-0.46) /
                      (2π * ε₀ * 10^(6 * 0.54))
        elseif identifier === :cigre2019
            σ = σ₀ + 4.7e-6 * σ₀^0.27 * f^0.54
            εr = 12 + 9.5e4 * σ₀^0.27 * f^(-0.46)
        elseif identifier === :datsios2019
            σ₄₂ = 1e4 / rho
            p = 0.537 * σ₄₂^0.16
            ε∞ = 1.24 * σ₄₂^0.415 * 3.5
            ε₃ = 4 * σ₄₂^0.463 * (2.9 * 3.5 - 3.8)
            fε = f < 3000 ? 3000 : f
            εr = ε∞ + (3000 / fε)^p * (ε₃ - ε∞)
            σ∞ = σ₄₂ * (1 + 0.65 / σ₄₂^0.57)
            σ = 1e-4 * (σ∞ * f * 1e-6 +
                (42.0 - 42.0 * (f - 42.0) * 1e-6) * (σ₄₂ / 42.0 - σ∞ * 1e-6))
        elseif identifier === :longmire1975
            fc = (125 * σ₀)^0.8312
            εr = 5 + sum(a / (1 + (f / (10^(n - 1) * fc))^2)
                for (n, a) in pairs(coefficients))
            σ = σ₀ + 2π * f * ε₀ * sum(a * (f / (10^(n - 1) * fc)) /
                (1 + (f / (10^(n - 1) * fc))^2) for (n, a) in pairs(coefficients))
        elseif identifier === :messier1985
            εr = 8 + sqrt(σ₀ * 8 / (π * f * ε₀))
            σ = σ₀ + sqrt(4π * f * σ₀ * ε₀ * 8)
        elseif identifier === :portela1999
            σ = σ₀ + 0.1e-6 * (2π * f)^0.72
            εr = 0.1e-6 * tan(π * 0.72 / 2) * (2π * f)^(-0.28) / ε₀
        elseif identifier === :scott1967
            u, x = log10(1000 / rho), log10(f)
            εr = 10^(5.491 + 0.946u - 1.097x + 0.069u^2 - 0.114u*x + 0.067x^2)
            σ̂ = 10.0^(0.028 + 1.098u - 0.068x + 0.036u^2 - 0.046u*x + 0.018x^2)
            return 1000.0 / σ̂, εr
        elseif identifier === :visacro1987
            εr = 2.34e6 * σ₀^0.535 * f^(-0.597)
            σ = σ₀ * (f / 100)^0.072
        elseif identifier === :visacro2012
            fe = f < 100 ? 100 : f
            εr = 1.3 + 7.6e3 * fe^(-0.4)
            σ = σ₀ + 1.2e-6 * σ₀^0.27 * (fe - 100)^0.65
        end
        return inv(σ), εr
    end

    causal = (:portela1999, :alipio2014, :longmire1975, :messier1985)
    identifiers = (:alipio2014, :cigre2019, :datsios2019, :longmire1975,
        :messier1985, :portela1999, :scott1967, :visacro1987, :visacro2012)
    for identifier in identifiers
        law = FD.Formula(identifier)
        max_rho, max_eps = 0.0, 0.0
        positive_rho, positive_eps = 0, 0
        symmetry_residual = 0.0
        for rho in (100.0, 1000.0, 10000.0)
            material = EP.EarthMaterial(rho, 10.0, 1.0)
            for f in 10.0 .^ range(0, 7; length = 61)
                ρ, εr = published(identifier, rho, f)
                evaluated = law(material, f)
                if identifier in causal
                    @test evaluated.rho ≈ ρ rtol = 1e-12
                    @test evaluated.eps_r ≈ εr rtol = 1e-12
                else
                    @test evaluated.rho == ρ
                    @test evaluated.eps_r == εr
                end
                max_rho = max(max_rho, abs(evaluated.rho / ρ - 1))
                max_eps = max(max_eps, abs(evaluated.eps_r / εr - 1))
            end

            s = 100.0 + 2π * 1e3 * im
            if identifier ∉ causal
                @test_throws ArgumentError law(material, 1e3; s)
                @test_throws Regex(string(identifier)) law(material, 1e3; s)
                continue
            end
            σ₀ = inv(rho)
            κ = if identifier === :portela1999
                σ₀ + (0.1e-6 / cos(π * 0.72 / 2)) * s^0.72
            elseif identifier === :alipio2014
                D = 1.26 * (1000 / rho)^(-0.73)
                σ₀ + s * ε₀ * 12 + (σ₀ * D / cos(π * 0.54 / 2)) *
                    (s / (2π * 1e6))^0.54
            elseif identifier === :longmire1975
                fc = (125 * σ₀)^0.8312
                σ₀ + s * ε₀ * 5 + sum(s * ε₀ * a /
                    (1 + s / (2π * 10^(n - 1) * fc)) for (n, a) in pairs(coefficients))
            elseif identifier === :messier1985
                σ₀ + s * ε₀ * 8 + 2 * sqrt(s * σ₀ * ε₀ * 8)
            end
            evaluated = law(material, 1e3; s)
            @test inv(evaluated.rho) + s * ε₀ * evaluated.eps_r ≈ κ rtol = 1e-12
            @test FD.constitutive(law, material, 1e3; s) == evaluated
            @test FD.constitutive(nothing, material, 1e3; s) === material
            conjugate = law(material, 1e3; s = conj(s))
            κ_conjugate = inv(conjugate.rho) + conj(s) * ε₀ * conjugate.eps_r
            κ_evaluated = inv(evaluated.rho) + s * ε₀ * evaluated.eps_r
            @test κ_conjugate ≈ conj(κ_evaluated) rtol = 1e-12
            symmetry_residual = max(symmetry_residual, abs(κ_conjugate - conj(κ_evaluated)))

            δ = 1.0
            points = (s + δ, s - δ, s + δ * im, s - δ * im)
            values = map(points) do point
                pair = law(material, imag(point) / 2π; s = point)
                inv(pair.rho) + point * ε₀ * pair.eps_r
            end
            dσ = (values[1] - values[2]) / (2δ)
            dω = (values[3] - values[4]) / (2δ)
            residual = abs(dσ + im * dω)
            println("Earth CR ", identifier, " rho0=", rho, " residual=", residual,
                " relative=", residual / abs(dσ))
            @test residual <= 1e-6 * abs(dσ)

            for f in (10.0, 1e3, 1e5), ratio in (0.1, 1.0, 10.0)
                point = 2π * f * (ratio + im)
                pair = law(material, f; s = point)
                positive_rho += pair.rho > 0
                positive_eps += pair.eps_r > 0
                @test pair.rho > 0
                @test pair.eps_r > 0
            end
        end
        println("Earth real-axis ", identifier, " max_relative_rho=", max_rho,
            " max_relative_eps=", max_eps)
        if identifier in causal
            println("Earth positivity summary ", identifier, " rho_positive=", positive_rho,
                "/27 eps_positive=", positive_eps, "/27")
            println("Earth conjugate symmetry ", identifier, " max_absolute_residual=", symmetry_residual)
            promoted = law(EP.EarthMaterial(100.0f0, 10.0f0, 1.0f0), 1e3;
                s = 100.0 + 2π * 1e3 * im)
            @test promoted == law(EP.EarthMaterial(100.0, 10.0, 1.0), 1e3;
                s = 100.0 + 2π * 1e3 * im)
        end
    end
end
