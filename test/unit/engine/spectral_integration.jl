@testitem "Engine / spectral integration algorithms and scalar contracts" tags=[:unit] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    options(method;
        controls...) = E.formulation_options(E.SpectralIntegral,
        (method = method, options = (; controls...)))
    using QuadGK
    a,b,h,y=1+im,2+im,0.5,1.0
    integral=E.SpectralIntegral(Val(:cosine),x->a*exp(-b*x),(height=h,separation=y),1.0)
    exact=a*(h+b)/((h+b)^2+y^2)
    q=1+.25im
    radial=E.SpectralIntegral(Val(:radial),u->a*exp(-b*u),(height=h,separation=y,q),1.0)
    # Independent hyperbolic parameterization. |cos(z)| <= exp(|Im(z)|).
    # For t>=T the integrand is bounded by |a|exp(-d*exp(t)),
    # d=(Re((h+b)q)-|Im(yq)|)/2; its integral <= exp(-d*exp(T))/(d*exp(T)).
    d=(real((h+b)*q)-abs(imag(y*q)))/2
    @test d>0
    T=log(80/d)
    radial_exact,estimate=quadgk(t->a*exp(-(h+b)*q*cosh(t))*cos(y*q*sinh(t)),
        0.0,T;rtol=1e-12,maxevals=10^6)
    uncertainty=estimate+abs(a)*exp(-d*exp(T))/(d*exp(T))
    for method in (:quad,:trapz,:cim)
        controls=method===:trapz ? (;max_refinements=14) : method===:cim ? (;samples=512,maxevals=10^6) : (;maxevals=10^6)
        o=options(method;controls...)
        @info "S4 resolved options" method o.options
        for (problem,expected,u) in ((integral,exact,0.0),(radial,radial_exact,uncertainty))
            actual=E.integrate(o.method,problem,o.options,nothing)
            for component in (real,imag)
                B=1e-6*abs(component(expected))
                @test u<=B/4
                @test abs(component(actual-expected))+u<=B
            end
        end
    end
    # The extracted pole is evaluated independently using its Laplace identity,
    # 1/(lambda+p)=integral(exp(-(lambda+p)t),t=0..Inf).
    residue,location=1+.2im,.5-.1im
    stop=80/real(location)
    pole_exact,pole_error=quadgk(t->residue*exp(-location*t)*(1+t)/((1+t)^2+.2^2),
        0.0,stop;rtol=1e-12,maxevals=10^6)
    pole_error+=abs(residue)*exp(-real(location)*stop)/(real(location)*(1+stop))
    for pole in (nothing,(;residue,location)),method in (:quad,:trapz,:cim)
        zero_residual=E.SpectralIntegral(Val(:cosine),lambda->zero(complex(lambda)),
            (height=1.0,separation=.2),1.0;pole)
        controls=method===:trapz ? (;max_refinements=14) : method===:cim ? (;samples=512,maxevals=10^6) : (;maxevals=10^6)
        o=options(method;controls...)
        actual=E.integrate(o.method,zero_residual,o.options,nothing)
        if pole===nothing
            @test iszero(actual)
        else
            # This is a quadrature control, with the approved S4 component
            # budget of 1e-6; the pole integral is not an algebraic identity.
            for component in (real,imag)
                budget=1e-6*abs(component(pole_exact))
                @test pole_error<=budget/4
                @test abs(component(actual-pole_exact))+pole_error<=budget
            end
        end
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
