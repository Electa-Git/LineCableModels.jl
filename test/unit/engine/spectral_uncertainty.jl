@testitem "Engine / spectral quadrature preserves correlated physical uncertainty" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    using Measurements
    const E=LineCableModels.Engine
    a=measurement(2.0, 0.02)
    rate=complex(3+0.2a, one(a))
    amplitude=complex(a, 0.3a)
    exact=amplitude/rate
    integral=E.SpectralIntegral(λ->amplitude*exp(-rate*λ))
    for method in (:quad,)
        controls=E.formulation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-9,)))
        value,
        error=E.integrate(controls.method, integral, controls.options; points = (0.0, 1.0))
        @test value isa Complex{Measurement{Float64}}
        @test error isa Float64
        @test E.numerical_magnitude(value-exact)<1e-9*E.numerical_magnitude(exact)
        for component in (real, imag)
            actual=Measurements.derivative(component(value), a)
            expected=Measurements.derivative(component(exact), a)
            @test actual≈expected rtol=1e-8
        end
    end
    # A zero nominal value must not erase the derivative during refinement.
    perturbation=measurement(0.0, 0.1)
    zero_mean=E.SpectralIntegral(λ->complex(perturbation)*exp(-(2+im)*λ))
    for method in (:quad,)
        controls=E.formulation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-8,)))
        actual,
        _=E.integrate(controls.method, zero_mean, controls.options; points = (0.0, 1.0))
        expected=complex(perturbation)/(2+im)
        @test E.numerical_magnitude(actual-expected)<1e-9
    end
    # Equal nominal media with independent uncertainties have no nominal
    # denominator pole. Metadata must not divide by their uncertain zero.
    sigma=[measurement(0.1, 0.001), measurement(0.1, 0.002)]
    s=complex(0.0, 2pi*1e4)
    epsilon=fill(8.8541878128e-12, 2)
    mu=fill(4pi*1e-7, 2)
    geometry=(horizontal = measurement.([0.0, 1.0]),
        height = measurement.([1.0, -1.0]), radius = measurement.([0.02, 0.03]))
    state=E.EarthImpedance._unified_state!(
        UnifiedFormulaFixtures.buffers(geometry).unified.current,
        (jω = s, Γ = zero(s), sigma, epsilon, mu),
        geometry)
    @test E.EarthImpedance.earth_denominator_pole(state)===nothing
    angle=E.EarthImpedance.earth_contour_angle(state, pi/6)
    @test all(isfinite,
        E.EarthImpedance._unified_points!(UnifiedFormulaFixtures.buffers(geometry).unified,
            state, 2.0, 1.0, 0.0, angle))
end
