@testitem "Engine / spectral integration algorithms and scalar contracts" tags=[:unit] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    options(method;
        controls...) = E.computation_options(E.SpectralIntegral,
        (method = method, options = (; controls...)))
    a, b, h, y = 2+im, 3+2im, 0.5, 1.2
    integral=E.SpectralIntegral(Val(:cosine), x->a*exp(-b*x), (height = h, separation = y), 1.0)
    exact=a*(h+b)/((h+b)^2+y^2)
    for method in (:quad, :trapz, :cim)
        o=options(method)
        @test E.integrate(o.method, integral, o.options, nothing) ≈
              exact rtol=2e-6
    end
    radial=E.SpectralIntegral(Val(:radial), u->a*exp(-b*u),
        (height = h, separation = y, q = 0.7+0.2im), 1.0)
    radial_exact=a*E.special_besselk(0, radial.weight.q*sqrt((h+b)^2+y^2))
    for method in (:quad, :trapz, :cim)
        o=options(method)
        @test E.integrate(o.method, radial, o.options, nothing) ≈
              radial_exact rtol=2e-6
    end
    # A zero residual and an exact extracted pole need no exponential images.
    for pole in (nothing, (residue = 1.0+0.2im, location = 0.5-0.1im))
        zero_residual=E.SpectralIntegral(Val(:cosine), λ->zero(complex(λ)),
            (height = 1.0, separation = 0.2), 1.0; pole)
        q=options(:quad);
        c=options(:cim)
        @test E.integrate(c.method, zero_residual, c.options, nothing) ≈
              E.integrate(q.method, zero_residual, q.options, nothing)
    end
    @test_throws ArgumentError options(:unknown)
    @test_throws ArgumentError options(:quad; max_terms = 10)
    @test_throws ArgumentError options(:trapz; rtol = -1)
    @test_throws ArgumentError options(:cim; samples = 3)
    for T in (Float32, Float64, BigFloat)
        x=E.SpectralIntegral(Val(:cosine), λ->complex(one(T))*exp(-λ),
            (height = one(T), separation = one(T)), one(T))
        for method in (:quad, :trapz)
            o=options(method; rtol = T===Float32 ? 1e-5 : 1e-6)
            result=E.integrate(o.method, x, o.options, nothing)
            @test result isa Complex{T}
            @test result ≈ T(2)/T(5) rtol=2e-5
        end
        if T===BigFloat
            o=options(:cim)
            @test_throws ArgumentError E.integrate(o.method, x, o.options, nothing)
        end
    end
    divergent=E.SpectralIntegral(
        Val(:cosine), λ->complex(1.0), (
            height = 0.0, separation = 0.0), 1.0)
    o=options(:cim; max_tail_refinements = 2)
    @test_throws ErrorException E.integrate(o.method, divergent, o.options, nothing)
end
