@testitem "Engine / spectral backends preserve correlated physical uncertainty" tags=[:unit] begin
    using Measurements
    const E=LineCableModels.Engine
    a=measurement(2.0, 0.02)
    rate=complex(3+0.2a, one(a))
    amplitude=complex(a, 0.3a)
    exact=amplitude/rate
    integral=E.SpectralIntegral(Val(:cosine), λ->amplitude*exp(-rate*λ),
        (height = 0.0, separation = 0.0), 1.0; features = E.SpectralFeatures([0.0, 1.0]))
    for method in (:quad, :trapz)
        controls=E.computation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-9,)))
        estimate=E.spectral_estimate(controls.method, integral, controls.options, nothing)
        @test estimate.value isa Complex{Measurement{Float64}}
        @test E.spectral_magnitude(estimate.value-exact)<1e-9*E.spectral_magnitude(exact)
        for component in (real, imag)
            actual=Measurements.derivative(component(estimate.value), a)
            expected=Measurements.derivative(component(exact), a)
            @test actual≈expected rtol=1e-8
        end
    end
    # A zero nominal value must not erase the derivative during refinement.
    perturbation=measurement(0.0, 0.1)
    zero_mean=E.SpectralIntegral(Val(:cosine), λ->complex(perturbation)*exp(-(2+im)*λ),
        (height = 0.0, separation = 0.0), 1.0; features = E.SpectralFeatures([0.0, 1.0]))
    for method in (:quad, :trapz)
        controls=E.computation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-8,)))
        actual=E.integrate(controls.method, zero_mean, controls.options, nothing)
        expected=complex(perturbation)/(2+im)
        @test E.spectral_magnitude(actual-expected)<1e-9
    end
    controls=E.computation_options(E.SpectralIntegral, (method = :cim,))
    @test_throws ArgumentError E.integrate(controls.method, integral, controls.options, nothing)
    uncertain_weight=E.SpectralIntegral(Val(:cosine), λ->complex(exp(-λ)),
        (height = measurement(1.0, 0.01), separation = 0.0), 1.0)
    @test_throws ArgumentError E.integrate(controls.method, uncertain_weight, controls.options, nothing)
    # Equal nominal media with independent uncertainties have no nominal
    # denominator pole. Metadata must not divide by their uncertain zero.
    sigma=[measurement(0.1, 0.001), measurement(0.1, 0.002)]
    s=complex(0.0, 2pi*1e4)
    epsilon=fill(8.8541878128e-12, 2)
    mu=fill(4pi*1e-7, 2)
    geometry=E.EarthReturnGeometry(
        measurement.([0.0, 1.0]), measurement.([1.0, -1.0]), measurement.([0.02, 0.03]))
    state=E.unified_earth_state(
        (jω = s, Γ = zero(s), sigma, epsilon, mu,
            gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon)),
        geometry)
    @test E.earth_denominator_pole(state)===nothing
    angle=E.earth_contour_angle(state, pi/6)
    @test all(isfinite, E.earth_spectral_features(state, 2.0, 1.0, 0.0, angle).points)
end
