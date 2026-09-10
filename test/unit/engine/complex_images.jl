@testitem "Engine / complex image certificates retain exact physical identities" tags=[:unit] begin
    const E=LineCableModels.Engine
    s=2pi*1e6im
    sigma=(0.0, 10.0)
    epsilon=(8.8541878128e-12, 8.8541878128e-12)
    mu=(4pi*1e-7, 4pi*1e-7)
    sh=sigma .+ s .* epsilon
    k2=s .* mu .* sh
    state=(s, Γ = zero(s), sh, mu, k2, k = E.outgoing_root.(k2))
    function integral(state, height, separation, logscale = 0.0)
        kernel=E.EarthRadialSpectrum{:Z, typeof(state), Float64}(state, logscale)
        features=E.earth_spectral_features(state, height, separation, 0.0, 0.3)
        E.SpectralIntegral(Val(:radial), kernel,
            (height, separation, q = state.k[2]), max(abs(state.k[2]), inv(height));
            angle = 0.3, features)
    end
    workspace=(cim = E.CIMWorkspace(), images = ComplexF64[], exponents = ComplexF64[])
    controls=E.computation_options(E.SpectralIntegral, (method = :cim,)).options
    quad=E.computation_options(E.SpectralIntegral, (method = :quad, options = (rtol = 1e-10,))).options
    distant=integral(state, 2.0, 2.0)
    first=E.spectral_estimate(Val(:cim), distant, controls, workspace)
    @test E.spectral_bounded(distant)
    @test first.samples>0
    @test first.value≈E.integrate(Val(:quad), distant, quad, nothing) rtol=1e-6
    fit_count=workspace.cim.statistics.fits[]
    certificates=workspace.cim.statistics.certifications[]
    repeated=@inferred E.spectral_estimate(Val(:cim), distant, controls, workspace)
    @test repeated.value==first.value
    @test repeated.evaluations==1
    @test repeated.samples==0
    manual=merge(controls,(samples=16,))
    @test E.spectral_estimate(Val(:cim),distant,manual,workspace).samples==0
    @test workspace.cim.statistics.fits[]==fit_count
    @test workspace.cim.statistics.certifications[]==certificates
    near=integral(state, 2.0, 0.0)
    reused=E.spectral_estimate(Val(:cim), near, controls, workspace)
    @test reused.value≈E.integrate(Val(:quad), near, quad, nothing) rtol=1e-6
    @test workspace.cim.statistics.fits[]==fit_count
    @test any(c->c.envelope, only(workspace.cim.fits).certificates)
    scaled=E.spectral_estimate(Val(:cim), integral(state, 2.0, 0.0, 0.7), controls, workspace)
    @test scaled.value≈exp(0.7)*reused.value rtol=2e-14
    @test workspace.cim.statistics.fits[]==fit_count
    @test scaled.error<=controls.rtol*abs(scaled.value)

    # Changed material values invalidate the image identity.
    changed_sh=(sh[1], 2sh[2])
    changed_k2=s .* mu .* changed_sh
    changed=merge(state, (sh = changed_sh, k2 = changed_k2, k = E.outgoing_root.(changed_k2)))
    @test E.cim_identity(integral(changed, 2.0, 0.0).kernel)!=E.cim_identity(near.kernel)
    @test E.cim_reuse_estimate(integral(changed, 2.0, 0.0), controls,
        workspace, Ref(0), ComplexF64)===nothing
    shifted=E.SpectralIntegral(Val(:radial), near.kernel,
        merge(near.weight, (q = 1.1near.weight.q,)), near.scale;
        angle = near.angle, features = near.features)
    @test E.cim_reuse_estimate(shifted, controls, workspace, Ref(0), ComplexF64)===nothing
    custom=E.SpectralIntegral(Val(:radial),near.kernel,near.weight,near.scale;
        angle=near.angle,features=E.SpectralFeatures(near.features.points;tail=x->exp(-x)))
    @test E.cim_reuse_estimate(custom,controls,workspace,Ref(0),ComplexF64)===nothing

    # Arbitrary closures have no declared immutable identity and are never
    # reused merely because the same callable object was passed again.
    amplitude=Ref(1.0)
    callback=E.SpectralIntegral(Val(:cosine), x->complex(amplitude[])*exp(-x),
        (height = 1.0, separation = 1.0), 1.0)
    one=E.integrate(Val(:cim), callback, controls, workspace)
    amplitude[]=2.0
    two=E.integrate(Val(:cim), callback, controls, workspace)
    @test two≈2one rtol=1e-6
    @test workspace.cim.statistics.fits[]==fit_count
end

@testitem "Engine / local trapezoid resolution retains oscillations and cancellation" tags=[:unit] begin
    const E=LineCableModels.Engine
    controls=E.computation_options(E.SpectralIntegral,
        (method = :trapz, options = (rtol = 1e-7,))).options
    make(y)=E.SpectralIntegral(Val(:cosine), x->complex(1.0),
        (height = 1.0, separation = y), 1.0)
    slow=E.spectral_trapz_points(make(0.0), controls, nothing)
    fast=E.spectral_trapz_points(make(50.0), controls, nothing)
    @test length(fast)>length(slow)
    @test last(fast)==1.0 # The tail remains part of the integral.
    estimate=E.spectral_estimate(Val(:trapz), make(50.0), controls, nothing)
    exact=complex(inv(1+50.0^2))
    @test estimate.value≈exact rtol=1e-7
    @test estimate.error<=controls.rtol*abs(estimate.value)
    @test abs(estimate.value-exact)<=max(estimate.error, 8eps(Float64)*abs(exact))
    @test_throws ArgumentError E.SpectralIntegral(Val(:besselcosine), x->complex(exp(-x)),
        (height = 1.0, separation = 0.0, radius = -1.0), 1.0)
end
