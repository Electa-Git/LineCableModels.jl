@testitem "Engine / EHEM / enforced layer selection and precision" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    frequencies=Float32[50, 500]
    model=EarthModel(80.0f0, 8.0f0, 1.0f0; thickness = 4.0f0)
    add!(model, EarthLayer(300.0f0, 12.0f0, 1.0f0))
    evaluated=EarthProperties.evaluate(CPEarth(), model, frequencies)

    bottom=EHEM.EnforceLayer()
    top=EHEM.EnforceLayer(layer = 2)
    middle=EHEM.EnforceLayer(layer = 3)
    @test description(bottom) == "Assume bottom layer"
    @test description(top) == "Assume top earth layer"
    @test description(middle) == "Assume layer 3"

    for (selector, source_index) in ((bottom, 3), (top, 2), (middle, 3))
        selected=selector(evaluated, model)
        rho, epsilon, mu=selected.rho, selected.epsilon, selected.mu
        @test size(rho) == (2, 2)
        @test eltype(rho) === Float32
        @test eltype(epsilon) === Float32
        @test eltype(mu) === Float32
        @test rho[1, :] == fill(model.layers[1].rho, 2)
        @test rho[2, :] == fill(model.layers[source_index].rho, 2)
        @test epsilon[2, :] == evaluated.epsilon[source_index, :]
        @test mu[2, :] == evaluated.mu[source_index, :]
    end

    @test_throws DomainError EHEM.EnforceLayer(layer = 0)
    @test_throws BoundsError EHEM.EnforceLayer(layer = 4)(evaluated, model)
end

@testitem "Engine / insulation formulations / analytical limits across precision" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics
] begin
    using TOML

    impedance_formulation=InsulationImpedance.Lossless()
    admittance_formulation=InsulationAdmittance.Lossless()
    @test description(impedance_formulation) == "Lossless insulation (ideal dielectric)"
    @test description(admittance_formulation) ==
          "Lossless insulation (ideal dielectric)"
    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "coaxial_capacitance.toml"
    ))["lossless_formulation"]

    for T in (Float32, Float64, BigFloat)
        setprecision(BigFloat, 128) do
            r_in=parse(T, "0.01")
            r_ex=parse(T, "0.02")
            relative_permeability=parse(T, "1.2")
            relative_permittivity=parse(T, "2.3")
            angular_frequency=T(2)*T(π)*T(50)
            s=Complex{T}(zero(T), angular_frequency)

            impedance=impedance_formulation(
                r_in,
                r_ex,
                relative_permeability,
                s
            )
            potential=admittance_formulation(
                r_in,
                r_ex,
                relative_permittivity,
                s,
                zero(T)
            )
            @test impedance isa Complex{T}
            @test potential isa Complex{T}
            @test TestNumerics.isapprox_scaled(
                impedance,
                Complex{T}(
                    zero(T),
                    T(parse(BigFloat, reference["series_impedance_imaginary"]))
                )
            )
            @test TestNumerics.isapprox_scaled(
                potential,
                Complex{T}(T(parse(BigFloat, reference["potential_coefficient"])))
            )
            @test iszero(impedance_formulation(zero(T), r_ex, one(T), s))
            @test iszero(impedance_formulation(r_ex, r_ex, one(T), s))
            @test iszero(admittance_formulation(zero(T), r_ex, one(T), s, zero(T)))
            @test iszero(admittance_formulation(r_ex, r_ex, one(T), s, zero(T)))
        end
    end
end

@testitem "Engine / internal impedance / passivity and solid-conductor limits" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    formulation=InternalImpedance.ScaledBessel()
    @test description(formulation) == "Scaled Bessel (Schelkunoff)"

    r_in=0.005
    r_ex=0.01
    rho=1.7241e-8
    relative_permeability=1.0
    s=ComplexF64(0.0, 2π*50.0)
    inner=formulation(:inner, r_in, r_ex, rho, relative_permeability, s)
    outer=formulation(:outer, r_in, r_ex, rho, relative_permeability, s)
    mutual=formulation(:mutual, r_in, r_ex, rho, relative_permeability, s)
    @test all(isfinite, (inner, outer, mutual))
    @test real(inner) > 0
    @test real(outer) > 0
    @test real(mutual) > 0
    @test imag(inner) >= 0
    @test imag(outer) >= 0

    @test iszero(formulation(:inner, 0.0, r_ex, rho, relative_permeability, s))
    @test iszero(formulation(:mutual, 0.0, r_ex, rho, relative_permeability, s))
    solid_outer=formulation(:outer, 0.0, r_ex, rho, relative_permeability, s)
    @test isfinite(solid_outer)
    @test real(solid_outer) > 0
    @test_throws ArgumentError formulation(
        :unsupported,
        r_in,
        r_ex,
        rho,
        relative_permeability,
        s
    )
end

@testitem "Engine / homogeneous earth / reciprocity and boundary validation" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    impedance_module=LineCableModels.Engine.EarthImpedance
    admittance_module=LineCableModels.Engine.EarthAdmittance
    epsilon0=8.8541878128e-12
    mu0=4π*1.0e-7
    rho=[Inf, 100.0]
    epsilon=[epsilon0, 10epsilon0]
    mu=[mu0, mu0]
    s=ComplexF64(0.0, 2π*50.0)

    impedance_cases=(
        (impedance_module.Papadopoulos(), [-1.0, -1.2]),
        (impedance_module.Pollaczek(), [-1.0, -1.2]),
        (impedance_module.Carson(), [10.0, 12.0]),
        (impedance_module.Papadopoulos(Γx = 1), [-1.0, -1.2])
    )
    @test description(first(impedance_cases)[1]) == "Papadopoulos"
    @test description(impedance_cases[2][1]) == "Pollaczek"
    @test description(impedance_cases[3][1]) == "Carson"
    for (formulation, heights) in impedance_cases
        mutual=formulation(:mutual, heights, 0.25, rho, epsilon, mu, s)
        reciprocal=formulation(:mutual, reverse(heights), 0.25, rho, epsilon, mu, s)
        self=formulation(:self, heights, 0.25, rho, epsilon, mu, s)
        @test isfinite(mutual)
        @test mutual ≈ reciprocal
        @test self == mutual
    end

    admittance_cases=(
        (admittance_module.Papadopoulos(), [-1.0, -1.2]),
        (admittance_module.Pollaczek(), [-1.0, -1.2]),
        (admittance_module.Images(), [10.0, 12.0]),
        (admittance_module.Papadopoulos(Γx = 1), [-1.0, -1.2])
    )
    @test description(first(admittance_cases)[1]) == "Papadopoulos"
    @test description(admittance_cases[2][1]) == "Pollaczek"
    @test description(admittance_cases[3][1]) == "Electrostatic images"
    for (formulation, heights) in admittance_cases
        mutual=formulation(:mutual, heights, 0.25, rho, epsilon, mu, s)
        reciprocal=formulation(:mutual, reverse(heights), 0.25, rho, epsilon, mu, s)
        @test isfinite(mutual)
        @test mutual ≈ reciprocal
        @test formulation(:self, heights, 0.25, rho, epsilon, mu, s) == mutual
    end

    underground_impedance=impedance_module.Papadopoulos()
    underground_admittance=admittance_module.Papadopoulos()
    for formulation in (underground_impedance, underground_admittance)
        @test_throws ArgumentError formulation(
            :unknown, [-1.0, -1.0], 0.2, rho, epsilon, mu, s)
        @test_throws ArgumentError formulation(:mutual, [-1.0], 0.2, rho, epsilon, mu, s)
        @test_throws ArgumentError formulation(
            :mutual, [0.0, -1.0], 0.2, rho, epsilon, mu, s)
        @test_throws ArgumentError formulation(
            :mutual, [1.0, 1.0], 0.2, rho, epsilon, mu, s)
    end
    @test_throws ArgumentError impedance_module._not(3)
    @test_throws ArgumentError admittance_module._not(3)
end
