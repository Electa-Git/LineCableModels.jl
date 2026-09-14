@testitem "Engine / declared spectral features and tails control adaptive sampling" tags=[:unit] begin
    const E=LineCableModels.Engine
    # Translate a narrow feature through unrelated initial grid positions, and
    # vary amplitude to exercise DE's amplitude-sensitive native estimator.
    for centre in (0.371, 9.001), amplitude in (1.0, 1e-8)

        width=centre<1 ? 0.003 : 0.0001
        points=[0.0; centre .+ width .* [-8, -2, -1, 0, 1, 2, 8]]
        feature=E.SpectralFeatures(points)
        kernel=λ->complex(amplitude)*exp(-((λ-centre)/width)^2)
        integral=E.SpectralIntegral(
            Val(:cosine), kernel, (height = 1.0, separation = 0.0), 1.0; features = feature)
        expected=amplitude*sqrt(pi)*width*exp(-centre+width^2/4)
        for method in (:quad, :trapz)
            controls=E.computation_options(E.SpectralIntegral, (
                method, options = (rtol = 1e-7,)))
            actual=E.integrate(controls.method, integral, controls.options, nothing)
            @test actual≈expected rtol=1e-6
        end
    end
    # The endpoint vanishes, while a later peak supplies the whole integral.
    # Its declared tail envelope is checked against independent integration.
    kernel=λ->complex(λ)*exp(-λ)
    features=E.SpectralFeatures([0.0, 1.0, 2.0]; tail = t->exp(-2t)*(t/2+1/4))
    integral=E.SpectralIntegral(
        Val(:cosine), kernel, (height = 1.0, separation = 0.0), 1.0; features)
    for method in (:quad, :trapz, :cim)
        controls=E.computation_options(E.SpectralIntegral, (
            method, options = (rtol = 1e-6,)))
        actual=E.integrate(controls.method, integral, controls.options, nothing)
        @test actual≈0.25 rtol=1e-5
    end
    invalid=E.SpectralIntegral(Val(:cosine), λ->complex(exp(-λ)),
        (height = 1.0, separation = 0.0), 1.0; features = E.SpectralFeatures([0.0, 1.0]; tail = t->0.0))
    controls=E.computation_options(E.SpectralIntegral, (method = :cim,))
    @test_throws ErrorException E.cim_true_tail(invalid, 1.0, 1.0, 1e-9, controls.options)
    narrow=E.SpectralIntegral(Val(:cosine), λ->complex(exp(-((λ-9)/1e-4)^2)),
        (height = 1.0, separation = 0.0), 1.0; features = E.SpectralFeatures(9 .+
                                                                             1e-4 .* [
            -8, -1, 0, 1, 8]))
    limited=E.computation_options(E.SpectralIntegral,
        (method = :cim,
            options = (samples = 16, max_terms = 8, max_refinements = 1)))
    # A small image budget may fail explicitly; it must not report zero merely
    # because the initial scale probes missed the declared narrow feature.
    try
        actual=E.integrate(Val(:cim), narrow, limited.options, nothing)
        @test actual≈sqrt(pi)*1e-4*exp(-9+1e-8/4) rtol=1e-5
    catch error
        @test error isa ErrorException
        @test occursin("CIM", sprint(showerror, error))||occursin(":cim", sprint(showerror, error))||
              occursin("sample budget",sprint(showerror,error))
    end
end

@testitem "Engine / radial validation retains small physical roots" tags=[:unit] begin
    const E=LineCableModels.Engine
    # The two roots differ by ten orders of magnitude. Reconstructing λ²
    # from ag²-κg² loses the entire air-root neighbourhood in Float64.
    state=(k2 = (complex(1e-20), complex(1.0)),
        k = (complex(1e-10), complex(1.0)), mu = (1.0, 3.0), sh = (1e-6im, 0.1+1e-5im))
    nodes=[0.0, 1e-10, 1e-9, 1e-6, 0.1]
    for kind in (:Z, :phi, :voltage)
        kernel=E.EarthRadialSpectrum{kind, typeof(state), Float64}(state, 0.0)
        count=Ref(0)
        radial=E.SpectralIntegral(Val(:radial), E.SpectralKernelCounter(kernel, count),
            (height = 2.0, separation = 1.0, q = state.k[2]), 1.0)
        original=E.SpectralIntegral(Val(:cosine),
            E.earth_spectrum(Val(kind), Val(2), Val(2), state, 1.0, 1.0, 0.0),
            (height = 2.0, separation = 1.0), 1.0)
        for λ in nodes
            @test radial(λ)≈original(λ) rtol=2e-14
        end
        @test count[]==length(nodes)
        _, values, _=E.cim_fit_samples(radial, nodes, 1.0, state.k[2])
        for (λ, value) in zip(nodes, values)
            ag=sqrt(λ^2+state.k2[2])
            @test value / ag * exp(-2ag) * cos(λ)≈original(λ) rtol=2e-14
        end
    end
end

@testitem "Engine / combined voltage paths equal the manuscript endpoint terms" tags=[:unit] begin
    const E=LineCableModels.Engine
    geometry=E.EarthReturnGeometry([0.0, 1.0], [1.2, -0.9], [0.01, 0.025])
    s=complex(0.0, 2pi*1e4)
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    u=E.unified_earth_state(
        (jω = s, Γ = 1e-4+2e-4im, sigma, epsilon, mu,
            gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon)),
        geometry)
    for P in (1, 2), Q in (1, 2), reference in (:deep, :interface)
        P==2 && reference===:deep && continue
        hp=abs(geometry.height[P])
        hq=abs(geometry.height[Q])
        radius=geometry.radius[P]
        padding=hq/2
        sq=u.scaling[Q]
        g=(; hp, hq, radius, padding, logscale = sq, i0minus = E.bessel_i0m1(u.x[P]))
        combined=E.EarthPathVoltageSpectrum{P, Q, reference, typeof(u), typeof(g)}(u, g)
        voltage=E.earth_spectrum(Val(:voltage), Val(P), Val(Q), u, hp, hq,
            u.scaling[P]+sq)
        endpoint=E.earth_spectrum(Val(:endpoint), Val(P), Val(Q), u, 0.0, hq, sq)
        surface=E.earth_spectrum(Val(:surface), Val(P), Val(Q), u, 0.0, hq, sq)
        for λ in (0.001, 0.1, 3.0, 40.0) .* exp(0.03im)
            j0=E.SpecialFunctions.besselj(0, radius*λ)
            original=u.A[P]*voltage(λ)*exp(-(hp+hq)*λ)
            P==1 && (original+=endpoint(λ)*j0*exp(-hq*λ))
            reference===:interface && (original-=surface(λ)*j0*exp(-hq*λ))
            actual=combined(λ)*exp(-(hq-padding)*λ)
            @test actual≈original rtol=2e-10 atol=1e-20
        end
        # At a_receiver=0, I0(κ_receiver*r)=J0(r*λ). The removable
        # quotient tends to -h_receiver*J0, including the spectral padding.
        roots=ntuple(m->m==P ? complex(-1.0) : u.k2[m], 2)
        branch=merge(u, (; k2 = roots))
        kernel=E.EarthPathVoltageSpectrum{P, Q, reference, typeof(branch), typeof(g)}(branch, g)
        a=map(k2->E.outgoing_root(1+k2), roots)
        other=a[3 - P]
        j0=E.SpecialFunctions.besselj(0, radius)
        expected=-hp*other*j0*exp(sq-hq*a[Q])/(u.sh[2]*a[1]+u.sh[1]*a[2])
        @test kernel(1.0)*exp(-(hq-padding))≈expected rtol=1e-12
    end
end

@testitem "Engine / prescribed longitudinal constants keep spectral contours on their branch" tags=[:unit] begin
    const E=LineCableModels.Engine
    geometry=E.EarthReturnGeometry([0.0, 1.0], [1.0, -1.0], [0.02, 0.03])
    s=complex(0.0, 2pi*1e4)
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    controls=E.computation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-9,))).options
    for Γ in (0.0+0.0im, 1e-4+2e-4im, 3e-4+1e-4im)
        state=E.unified_earth_state(
            (jω = s, Γ, sigma, epsilon, mu,
                gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon)),
            geometry)
        angle=E.earth_contour_angle(state, pi/6)
        @test 0<angle<=pi/6
        iszero(Γ) || @test angle<pi/6
        for P in (1, 2), Q in (1, 2), kind in (:Z, :phi, :voltage)
            # At Γ=0 the separate air-voltage term has a branch-point
            # singularity; production combines its cancelling endpoints.
            # Here the added zero-Γ check concerns the scalar coordinate.
            iszero(Γ)&&kind!==:phi && continue
            kernel=E.earth_spectrum(Val(kind), Val(P), Val(Q), state, 1.0, 1.0, 0.0)
            integral=E.SpectralIntegral(
                Val(:cosine), kernel, (height = 2.0, separation = 1.0), 1.0;
                features = E.earth_spectral_features(state, 2.0, 1.0, 0.0, 0.0))
            real_axis=E.integrate(Val(:quad), integral, controls, nothing)
            rotated=E.earth_spectral_term(
                Val(kind), Val(P), Val(Q), state, 1.0, 1.0, 1.0, 0.0, 0.0,
                Val(:quad), controls, nothing)
            @test rotated.value≈real_axis rtol=1e-7
            if kind===:phi
                cim=E.computation_options(E.SpectralIntegral,
                    (method = :cim, options = (rtol = 1e-7,))).options
                images=E.earth_spectral_term(Val(:phi), Val(P), Val(Q), state,
                    1.0, 1.0, 1.0, 0.0, 0.0, Val(:cim), cim, nothing)
                @test images.value≈real_axis rtol=1e-6
                @test images.error<=1e-6*abs(real_axis)
            end
        end
    end
end
