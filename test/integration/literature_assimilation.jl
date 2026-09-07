@testitem "Engine / literature assimilation / isolated internal impedances" begin
    using LineCableModels
    using SpecialFunctions: besseli
    II = LineCableModels.Engine.InternalImpedance

    for T in (Float32, Float64, BigFloat)
        r, rho, mur = T(1) / 100, T(17) / 10^9, one(T)
        μ = T(4) * (one(T) * π) / T(10)^7
        s = complex(zero(T), T(100) * (one(T) * π))
        for id in (:Knight2016, :BrandaoFaria2011)
            @test id in II.formulas()
            f = II.Formula(id)
            z = f(zero(T), r, rho, mur, s)(Val(:outer))
            @test z isa Complex{T}
            @test real(z) > 0 && imag(z) > 0
            @test f(zero(T), r, rho, mur, -s)(Val(:outer)) ≈ conj(z)
            @test iszero(f(zero(T), r, rho, mur, s)(Val(:inner)))
            @test iszero(f(zero(T), r, rho, mur, s)(Val(:mutual)))
        end
        @test II.Formula(:Knight2016)(zero(T), r, rho, mur, zero(s))(Val(:outer)) ≈
              rho / ((one(T) * π) * r^2)
        @test_throws DomainError II.Formula(:Knight2016)(r / 2, r, rho, mur, s)
        @test_throws DomainError II.Formula(:Knight2016)(zero(T), r, rho, mur, s + one(T))
        for p in T.((-3, 0, 2)), ratio in T.((0.1, 0.8))
            a = ratio * r
            f = II.Formula(:BrandaoFaria2011; exponent = p)
            z = f(a, r, rho, mur, s)(Val(:outer))
            root = sqrt((p / 2)^2 + s * μ * r^2 / rho)
            m1, m2 = p / 2 + root, p / 2 - root
            reference = rho / (2 * (one(T) * π) * r^2) *
                        (m2 * ratio^m2 - m1 * ratio^m1) /
                        (ratio^m1 - ratio^m2)
            @test z ≈ reference rtol = max(T(1e-12), 32eps(T))
            @test f(zero(T), r, rho, mur, s)(Val(:outer)) ≈ s * μ / (2 * (one(T) * π) * m1)
            dc = iszero(p) ? rho / (2 * (one(T) * π) * r^2 * log(inv(ratio))) :
                 rho * p / (2 * (one(T) * π) * r^2 * (ratio^(-p) - 1))
            @test f(a, r, rho, mur, zero(s))(Val(:outer)) ≈ dc
            @test_throws ArgumentError f(a, r, rho, mur, s)(Val(:inner))
            @test_throws ArgumentError f(a, r, rho, mur, s)(Val(:mutual))
        end
        f = II.Formula(:BrandaoFaria2011; exponent = 2)
        @test isfinite(f(r / 2, r, rho, mur, s * T(10)^12)(Val(:outer)))
        @test isfinite(f(r / 2, r, rho, mur, s / T(10)^12)(Val(:outer)))
    end
    @test_throws DomainError II.Formula(:BrandaoFaria2011; exponent = Inf)
    @test_throws ArgumentError II.Formula(:BrandaoFaria2011; unknown = 1)

    # Compare the fitted resistance and inductance separately with the
    # published cylindrical Bessel solution, across both asymptotes.
    knight = II.Formula(:Knight2016)
    schelkunoff = II.Formula(:Schelkunoff1934)
    for frequency in 10.0 .^ range(-8, 12; length = 200)
        s = complex(0.0, 2π * frequency)
        fitted = knight(0.0, 0.01, 1.7e-8, 1.0, s)(Val(:outer))
        exact = schelkunoff(0.0, 0.01, 1.7e-8, 1.0, s)(Val(:outer))
        @test abs(real(fitted) / real(exact) - 1) < 0.0009
        @test abs(imag(fitted) / imag(exact) - 1) < 0.00016
    end

    # An isolated tube needs no transfer term in the coaxial assembly.
    conductor = Material(kind = :conductor, rho = 1.7e-8, eps_r = 1.0,
        mu_r = 1.0, T0 = 20.0, alpha = 0.0)
    for shape in (Disk(0.01), Annulus(0.005, 0.01))
        design = build(CableDesign, "isolated", Stack(
            Group(:core, Region(:metal, shape, conductor))))
        engine = LineCableModels.Engine
        data = engine.LocalCableData(engine.flatten(LineCableModelsCoaxial(), design))
        method = II.Formula(:BrandaoFaria2011; exponent = -1)
        methods = (internal_impedance = method,
            insulation_impedance = engine.InsulationImpedance.Formula(:Ametani1980))
        destination = zeros(ComplexF64, 1, 1)
        engine.cable_impedance!(destination, data, [1.7e-8], methods, complex(0.0, 100π))
        @test destination[1, 1] ≈ method(data.r_in[1], data.r_ext[1],
            1.7e-8, 1.0, complex(0.0, 100π))(Val(:outer))
    end
end

@testitem "Engine / literature assimilation / exact source frontmatter" begin
    using Base.Docs: meta
    using LineCableModels
    root = pkgdir(LineCableModels)
    owner = LineCableModels.Engine.InternalImpedance
    records = (
        :Knight2016 => "internal-impedance/2016/solid-round-continuous-fitted-approximation/Knight2016.md",
        :BrandaoFaria2011 => "internal-impedance/2011/inhomogeneous-euler-cauchy-tubular-conductor/BrandaoFaria2011.md"
    )
    for (id, path) in records
        source = read(joinpath(root, "docs", "theory", path), String)
        source_frontmatter = match(r"(?s)(## Identification and source\n.*?)(?=\n\n\*\*Description\.\*\*)", source)[1]
        matches = String[]
        for (binding, multidoc) in meta(owner)
            binding.var === :description || continue
            for (typesig, docstring) in multidoc.docs
                first(Base.unwrap_unionall(only(typesig.parameters)).parameters) === id || continue
                push!(matches, join(filter(x -> x isa String, collect(docstring.text))))
            end
        end
        @test length(matches) == 1
        @test occursin(source_frontmatter, only(matches))
    end
end
