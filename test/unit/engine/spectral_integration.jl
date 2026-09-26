@testitem "Engine / callable-only QuadGK integration and scalar behavior" tags=[:unit] begin
    const E=LineCableModels.Engine
    options(;
        controls...)=E.formulation_options(E.SpectralIntegral,
        (method = :quad, options = (; controls...)))
    calls=Ref(0)
    f=x->begin
        calls[]+=1
        (1+im)*exp(-(1.5+im)*x)*cos(x)
    end
    integral=E.SpectralIntegral(f)
    @test calls[]==0
    @test fieldnames(typeof(integral))==(:f,)
    @test integral.f===f
    settings=options()
    @test settings.method===Val(:quad)
    actual, error=E.integrate(settings.method, integral, settings.options)
    exact=(1+im)*(1.5+im)/((1.5+im)^2+1)
    @test actual≈exact rtol=1e-8
    @test error>=0 && isfinite(error)
    @test_throws ArgumentError E.formulation_options(E.SpectralIntegral, (method = :unknown,))
    for controls in ((max_terms = 10,), (rtol = -1,), (maxevals = 0,), (maxevals = true,),
        (rtol = 0, atol = 0), (samples = 16,))
        @test_throws ArgumentError options(; controls...)
    end
    for method in (:trapz, :cim)
        @test_throws ArgumentError E.formulation_options(E.SpectralIntegral, (; method))
    end
    for T in (Float32, Float64, BigFloat)
        scalar_settings=options(rtol = T===Float32 ? 1e-5 : 1e-9)
        scalar_integral=E.SpectralIntegral(x->complex(exp(-2x)*cos(x)))
        scalar_value,
        scalar_error=E.integrate(
            scalar_settings.method, scalar_integral, scalar_settings.options;
            coordinate_type = T, points = T[1])
        @test scalar_value isa Complex{T}
        @test scalar_error isa T
        @test scalar_value≈T(2)/T(5) rtol=2e-5
    end
    @test_throws DomainError E.integrate(Val(:quad), E.SpectralIntegral(x->NaN), settings.options)
end

@testitem "Engine / real-only callable struct and complete transformed integrands" tags=[:unit] begin
    const E=LineCableModels.Engine
    struct RealKernel
        rate::Float64
    end
    (k::RealKernel)(x::Real)=exp(-k.rate*x)
    controls=(rtol = 1e-10, atol = 0.0, maxevals = 100000)
    value, error=E.integrate(Val(:quad), E.SpectralIntegral(RealKernel(2.0)), controls)
    @test value isa Float64
    @test value≈0.5 rtol=1e-10
    @test error<=controls.rtol*abs(value)
    for c in (0.1, 2.0, 3cis(0.2))
        f=t->c*exp(-2c*t)
        transformed,
        estimate=E.integrate(Val(:quad), E.SpectralIntegral(f), controls; points = (1.0,))
        @test transformed≈0.5 rtol=1e-10
        @test estimate<=controls.rtol*abs(transformed)
    end
end

@testitem "Engine / QuadGK reports an unmet target without retry or rejection" tags=[:unit] begin
    using QuadGK: quadgk
    const E=LineCableModels.Engine
    calls=Ref(0)
    f=x->begin
        calls[]+=1
        exp(-x)*cos(40x)
    end
    controls=(rtol = 1e-13, atol = 0.0, maxevals = 15)
    actual=@test_logs (:warn, r"QuadGK returned an estimated error") E.integrate(
        Val(:quad), E.SpectralIntegral(f), controls; points = (1.0,))
    work=calls[]
    calls[]=0
    mapped=t->begin
        den=inv(1-t)
        f(t*den)*den^2
    end
    direct=quadgk(mapped, 0.0, 0.5, 1.0; controls..., norm = E.numerical_magnitude)
    @test work==calls[]
    @test actual==direct
    @test isfinite(first(actual))
    @test last(actual)>controls.rtol*abs(first(actual))
end

@testitem "Engine / complete contour integrands preserve half-line endpoints" tags=[:unit] begin
    const E=LineCableModels.Engine
    controls=(rtol = 1e-10, atol = 0.0, maxevals = 100000)
    # The rational map has a nonzero limiting value at the infinite endpoint.
    for angle in (0.0, 0.2), points in ((), Float64[])

        c=cis(angle)
        integral=E.SpectralIntegral(t->c/(1+c*t)^2)
        value, error=E.integrate(Val(:quad), integral, controls; points)
        @test value≈1 rtol=1e-10
        @test abs(value-1)<=error+4eps(Float64)
    end
end

@testitem "Engine / quadrature normalizes numerical points and reuses storage" tags=[:unit] begin
    const E=LineCableModels.Engine
    for T in (Float32, Float64, BigFloat)
        controls=(rtol = T===Float32 ? 1e-5 : 1e-9, atol = 0.0, maxevals = 100000)
        points=T[2, 0, 1, 2, 1]
        retained=copy(points)
        integral=E.SpectralIntegral(x->complex(exp(-x)))
        workspace=E.integration_workspace(T, Complex{T})
        expected=E.integrate(Val(:quad), integral, controls; points, coordinate_type = T)
        for _ in 1:2
            value, error=E.integrate(Val(:quad), integral, controls, workspace; points)
            @test value isa Complex{T}
            @test value≈one(T) rtol=controls.rtol
            @test value≈first(expected) rtol=controls.rtol
            @test error<=controls.rtol*abs(value)
            @test workspace.points==T[0, 1, 2]
        end
        @test points==retained
        for invalid in (T[-1], T[Inf], T[NaN], [1im])
            @test_throws DomainError E.integrate(
                Val(:quad), integral, controls, workspace; points = invalid)
        end
    end
end

@testitem "Engine / complete oscillatory and Bessel integrands" tags=[:unit] begin
    using SpecialFunctions: besselj
    const E=LineCableModels.Engine
    controls=(rtol = 1e-7, atol = 0.0, maxevals = 1000000)
    oscillatory=E.SpectralIntegral(x->complex(exp(-x)*cos(50x)))
    value, error=E.integrate(Val(:quad), oscillatory, controls; points = (1.0,))
    exact=complex(inv(1+50.0^2))
    @test value≈exact rtol=1e-7
    @test error<=controls.rtol*abs(value)
    @test abs(value-exact)<=max(error, 8eps(Float64)*abs(exact))
    bessel=E.SpectralIntegral(x->complex(exp(-x)*besselj(0, 0.3x)))
    @test first(E.integrate(Val(:quad), bessel, controls))≈inv(sqrt(1+0.3^2)) rtol=1e-7
    # The zero-radius identity belongs to the supplied expression, not to a
    # special numerical-executor branch.
    f=x->exp(-x)*cos(0.2x)/complex(x+0.01)
    zero_radius=E.SpectralIntegral(x->f(x)*besselj(0, zero(x)))
    for x in (0.0, 0.01, 1.0)
        @test zero_radius.f(x)==f(x)
    end
    @test first(E.integrate(Val(:quad), zero_radius, controls; points = (0.01, 1.0))) ≈
          first(E.integrate(Val(:quad), E.SpectralIntegral(f), controls; points = (
        0.01, 1.0))) rtol=1e-7
end

@testitem "Engine / explicit numerical points resolve displaced narrow peaks" tags=[:unit] begin
    const E=LineCableModels.Engine
    controls=(rtol = 1e-7, atol = 0.0, maxevals = 1000000)
    for centre in (0.371, 9.001), amplitude in (1.0, 1e-8)

        width=centre<1 ? 0.003 : 0.0001
        points=[0.0; centre .+ width .* [-8, -2, -1, 0, 1, 2, 8]; 1.0]
        integral=E.SpectralIntegral(x->complex(amplitude)*exp(-((x-centre)/width)^2-x))
        expected=amplitude*sqrt(pi)*width*exp(-centre+width^2/4)
        actual, _=E.integrate(Val(:quad), integral, controls; points)
        @test actual≈expected rtol=1e-6
    end
    # A vanishing endpoint is not a license to discard a nonzero integral.
    value,
    _=E.integrate(Val(:quad), E.SpectralIntegral(x->complex(x)*exp(-2x)),
        controls; points = (0.0, 1.0, 2.0))
    @test value≈0.25 rtol=1e-5
end
