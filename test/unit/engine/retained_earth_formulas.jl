@testitem "Engine / retained earth equations preserve numerical baselines across integration methods" tags=[:unit] begin
    using TOML
    const E=LineCableModels.Engine
    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels), "test/fixtures/reference/retained_earth_formulas.toml"))
    rho=[Inf, 100.0]
    epsilon=8.8541878128e-12 .* [1, 10]
    mu=fill(4pi*1e-7, 2)
    for row in reference["cases"]
        owner=getproperty(E, Symbol(row["family"]))
        id=Symbol(row["id"])
        placement=Symbol(row["placement"])
        self=row["kind"]=="self"
        heights=placement===:overhead ? (1.0, 2.0) :
                placement===:underground ? (-1.0, -2.0) : (1.0, -2.0)
        self && (heights=(heights[1], heights[1]))
        layers=placement===:overhead ? (1, 1) : placement===:underground ? (2, 2) : (1, 2)
        pair=E.EarthPair(1, self ? 1 : 2, heights, self ? 0.0 : 0.75,
            layers; radius = self ? 0.01 : nothing)
        s=complex(0.0, 2pi*row["frequency"])
        expected=complex(row["real"], row["imag"])
        if id === :Pollaczek1926 && placement !== :underground
            @test_throws ArgumentError owner.Formula(id)(rho, epsilon, mu, s, pair)
            continue
        end
        binding = validate(owner.Formula(id), pair)
        methods = haskey(binding.options, :integration) ? (:quad, :trapz, :cim) : (:direct,)
        for method in methods
            @testset "$(row["family"]) $id $placement $(row["kind"]) $(row["frequency"]) $method" begin
                options = method === :direct ? (;) :
                          (integration = (method = method, options = (;)),)
                actual=owner.Formula(id; options)(rho, epsilon, mu, s, pair)()
                @test actual≈expected rtol=(method===:quad ? 2e-8 : 3e-6) atol=0
                # Compare the small dissipative component separately.
                @test real(actual)≈real(expected) rtol=(method===:quad ? 2e-7 : 3e-5) atol=1e-14
            end
        end
    end
end

@testitem "Engine / overhead potential pole extraction agrees with the unsplit real-axis equation" tags=[:unit] begin
    using QuadGK
    const E=LineCableModels.Engine
    # This reference uses the original spectral expression, independent of
    # SpectralIntegral's analytic pole extraction and contour implementation.
    epsilon=8.8541878128e-12 .* [1, 10]
    mu=fill(4pi*1e-7, 2)
    rho=[Inf, 100.0]
    for freq in (50.0, 1e5), separation in (0.75, 10.0)

        s=complex(0.0, 2pi*freq)
        g0=s^2*mu[1]*epsilon[1]
        g1=s*mu[2]*(inv(rho[2])+s*epsilon[2])
        r=g1/g0
        a=sqrt(g1-g0)
        H=3.0
        kernel(λ) = exp(-H*λ)*cos(separation*λ)/(r*λ+sqrt(λ^2+g1-g0))
        original=first(quadgk(kernel, 0.0, abs(a/r), abs(a), 1.0, Inf; rtol = 1e-10))
        pair=E.EarthPair(1, 2, (1.0, 2.0), separation, (1, 1))
        direct=log(hypot(separation, H)/hypot(separation, 1.0))
        for method in (:quad, :trapz, :cim)
            o=(integration = (method = method, options = (;)),)
            coefficient=E.EarthAdmittance.Formula(:Wise1948; options = o)(
                rho, epsilon, mu, s, pair)()
            recovered=(coefficient*(2pi*epsilon[1])-direct)/2
            @test recovered≈original rtol=2e-6 atol=1e-14
        end
    end
end
