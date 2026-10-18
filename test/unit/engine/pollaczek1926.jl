@testitem "Engine / Pollaczek underground integral / independent values and materials" tags=[:unit, :engine] begin
    const E = LineCableModels.Engine
    # Independent mpmath 65-digit quadrature of Saad (1996), Eqs. (1)-(4),
    # with the even full-line J reduced to twice the half-line cosine integral.
    fixtures = (
        ((-1, -1), 0, 0.01, 1, 0,
            "0.0000494646829716793566621026239403048913564914585774273984247",
            "0.0007188215777546735924839215997720618297443189025802483192512"),
        ((-1, -2), 0.75, nothing, 2.5, 25,
            "0.00007786307992028348054111885965150945848404391152360497993178",
            "0.0003826434452327189821792578688421453965267331776601825445374"))
    for T in (Float32, Float64, BigFloat)
        μ0, ε0 = 4T(pi)*T(10)^(-7), parse(T, "8.8541878128e-12")
        selected = E.EarthImpedance.Formula(:pollaczek1926;
            options=(integration=(options=(rtol=T === Float32 ? 1e-5 : 1e-12,),),))
        for (heights, x, radius, μr, sigma, re, imz) in fixtures
            pair = E.EarthPair(1, radius === nothing ? 2 : 1, T.(heights), T(x), (2, 2);
                radius=radius === nothing ? nothing : T(radius))
            s = complex(T(sigma), 100T(pi))
            rho, epsilon, mu = T[Inf, 100], [ε0, 10ε0], [μ0, T(μr)*μ0]
            value = @inferred selected(rho, epsilon, mu, s, pair)
            @test value isa Complex{T}
            @test isfinite(value)
            T === Float32 || (@test value ≈ complex(parse(T, re), parse(T, imz)) rtol=5e-13)
            @test selected(rho, 7epsilon, mu, s, pair) == value
            # Only the ratio mu_earth/rho_earth enters m. The prefactor uses mu_0.
            @test selected(T[Inf, 200], epsilon, [μ0, 2mu[2]], s, pair) == value
        end
    end
    # Independent 65-digit quadrature also resolves the largest Saad discrepancy:
    # at 100 ohm*m and 1 MHz, the mutual approximation differs by 2.160912%.
    selected = E.EarthImpedance.Formula(:pollaczek1926;
        options=(integration=(options=(rtol=1e-12,),),))
    pair = E.EarthPair(1, 2, (-1.0, -1.0), 1.0, (2, 2))
    value = selected([Inf, 100.0], [8.8541878128e-12, 8.8541878128e-11],
        [4pi*1e-7, 4pi*1e-7], complex(0.0, 2pi*1e6), pair)
    @test value ≈ complex(1.103498320106650833647385273886917874708014798129,
        2.095070752403021477192139482804206890454136563618) rtol=5e-13

end

@testitem "Engine / Pollaczek underground integral / convergence, asymptote and complex frequency" tags=[:unit, :engine] begin
    using Base.MathConstants: eulergamma
    const E = LineCableModels.Engine
    μ0, ε0 = 4pi*1e-7, 8.8541878128e-12
    selected = E.EarthImpedance.Formula(:pollaczek1926;
        options=(integration=(options=(rtol=1e-7,),),))
    tight = E.EarthImpedance.Formula(:pollaczek1926;
        options=(integration=(options=(rtol=1e-12,),),))
    pairs = (E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius=0.08),
        E.EarthPair(1, 2, (-1.0, -1.5), 1.0, (2, 2)))
    maxima = zeros(2)
    for pair in pairs, rho in (100.0, 1000.0, 10000.0)
        evaluate(formula, s) = formula([Inf, rho], [ε0, 10ε0], [μ0, μ0], s, pair)
        for f in (1e-3, 1.0, 1e3, 1e6, 5e6), sigma in (0.0, 80.0)
            s = complex(sigma, 2pi*f)
            value = @inferred evaluate(tight, s)
            difference = abs(evaluate(selected, s)-value)/abs(value)
            maxima[1] = max(maxima[1], difference)
            @test difference < 1e-7
            @test isfinite(value)
            @test evaluate(tight, conj(s)) ≈ conj(value) rtol=1e-12
            @test abs(evaluate(tight, s*(1+1e-7))-value)/abs(value) < 1e-6
        end
        # Derive the buried Carson-type limit from Saad (1996), Eqs. (1)-(4):
        # lambda=m*sinh(t) gives J=K0(m*D)+integral(exp(-m*ell*cosh(t))*
        # cos(m*x*sinh(t))*exp(-2t), t=0..Inf). The last integral has limit 1/2.
        # This gives Z/(s*mu0/(2pi))=-log(m*d/2)-EulerGamma+1/2+o(1), independent
        # of burial depth as s -> 0. For s=j*w, Re(Z) -> w*mu0/8.
        residuals = Float64[]
        d = hypot(something(pair.radius, pair.separation), pair.heights[1]-pair.heights[2])
        for f in (1e-8, 1e-12)
            s = complex(0.0, 2pi*f)
            m = sqrt(s*μ0/rho)
            value = evaluate(tight, s)
            residual = abs(value/(s*μ0/(2pi)) + log(m*d/2) + eulergamma - 0.5)
            push!(residuals, residual)
            @test residual < 1e-6
            @test real(value) ≈ imag(s)*μ0/8 rtol=1e-6
        end
        @test residuals[2] < residuals[1]
        maxima[2] = max(maxima[2], last(residuals))
    end
    forward, reverse = pairs[2], E.EarthPair(2, 1, (-1.5, -1.0), 1.0, (2, 2))
    for s in (80.0+0.0im, 80.0+2pi*im, 80.0+2pi*1e6*im)
        @test tight([Inf, 100.0], [ε0, 10ε0], [μ0, μ0], s, forward) ==
              tight([Inf, 100.0], [ε0, 10ε0], [μ0, μ0], s, reverse)
    end
    real_axis = tight([Inf, 100.0], [ε0, 10ε0], [μ0, μ0], 80.0+0.0im, forward)
    @test isfinite(real_axis)
    @test iszero(imag(real_axis))
    @info "Pollaczek quadrature and logarithmic limit" maximum_relative_change=maxima[1] maximum_low_frequency_bracket_error=maxima[2]
end

@testitem "Engine / Pollaczek underground integral / Saad benchmark comparison" tags=[:unit, :engine] setup=[TestFixtures] begin
    const E = LineCableModels.Engine
    problem = TestFixtures.known_case_problem(:cable_320kv_armoured_dc_bipole; frequencies=[1.0])
    radius = outer_radius(first(problem.system.designs))
    pairs = (E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius),
        E.EarthPair(1, 2, (-1.0, -1.0), 1.0, (2, 2)))
    selected = E.EarthImpedance.Formula(:pollaczek1926)
    saad = E.EarthImpedance.Formula(:saad1996)
    μ0, ε0 = 4pi*1e-7, 8.8541878128e-12
    for rho in (100.0, 1000.0, 10000.0)
        differences = map(pairs) do pair
            maximum(10.0 .^ range(0, 6; length=241)) do f
                arguments = ([Inf, rho], [ε0, 10ε0], [μ0, μ0], complex(0.0, 2pi*f), pair)
                exact = selected(arguments...)
                abs(saad(arguments...)-exact)/abs(exact)
            end
        end
        # The self differences stay below Saad's typical ~1.5%. The mutual
        # difference at 100 ohm*m reaches 2.16% at 1 MHz. The independent
        # integral fixture confirms this value, within the first kernel's ~3%.
        @test differences[1] < 0.015
        @test differences[2] < 0.03
        @info "Pollaczek versus Saad, 1 Hz to 1 MHz" rho radius self_percent=100differences[1] mutual_percent=100differences[2]
    end
end
