@testitem "Engine / analytic earth tails cover ordered media and voltage paths" tags=[:unit] begin
    const E=LineCableModels.Engine
    geometry=E.EarthReturnGeometry([0.0,1.0],[1.2,-0.9],[0.01,0.025])
    for frequency in (1.0,1e4,1e6), Γ in (0.0im,1e-4+2e-4im)
        s=2pi*frequency*im
        sigma=[0.0,0.1]; epsilon=8.8541878128e-12.*[1.0,8.0]; mu=4pi*1e-7.*[1.0,3.0]
        state=E.unified_earth_state((jω=s,Γ,sigma,epsilon,mu,
            gamma_medium_squared=s.*mu.*(sigma.+s.*epsilon)),geometry)
        angle=E.earth_contour_angle(state,pi/8)
        rotation=cis(angle)
        for P in 1:2, Q in 1:2
            hp,hq=1.2,0.9
            features=E.earth_spectral_features(state,hp+hq,1.0,0.02,angle)
            for kind in (:Z,:phi,:voltage,:endpoint,:surface,:finite_direct,:finite_image)
                kernel=E.earth_spectrum(Val(kind),Val(P),Val(Q),state,hp,hq,0.3)
                integral=E.SpectralIntegral(Val(:besselcosine),kernel,
                    (height=hp+hq,separation=1.0,radius=0.02),1.0;angle,features)
                @test E.spectral_bounded(integral)
                limit=max(20.0,4maximum(abs,state.k))
                bound=E.spectral_tail(integral)(limit)
                remainder,error=E.quadgk(t->abs(integral(rotation*(limit+t/(1-t))))/(1-t)^2,
                    0.0,1.0;rtol=1e-7,atol=max(bound*1e-8,floatmin(Float64)))
                @test remainder<=bound+error
                @test E.spectral_tail(integral)(2limit)<=bound
            end
            for reference in (:deep,:interface)
                padding=hq/2
                g=(hp,hq,radius=0.025,padding,logscale=0.3,i0minus=E.bessel_i0m1(state.k[P]*0.025))
                kernel=E.EarthPathVoltageSpectrum{P,Q,reference,typeof(state),typeof(g)}(state,g)
                path_features=E.earth_spectral_features(state,hq-padding,1.0,0.025,angle)
                integral=E.SpectralIntegral(Val(:cosine),kernel,
                    (height=hq-padding,separation=1.0),1.0;angle,features=path_features)
                limit=max(30.0,4maximum(abs,state.k))
                bound=E.spectral_tail(integral)(limit)
                remainder,error=E.quadgk(t->abs(integral(rotation*(limit+t/(1-t))))/(1-t)^2,
                    0.0,1.0;rtol=1e-7,atol=max(bound*1e-8,floatmin(Float64)))
                @test isfinite(bound) && remainder<=bound+error
            end
        end
    end
end

@testitem "Engine / radial and transformed tails bound the physical contour" tags=[:unit] begin
    const E=LineCableModels.Engine
    for contrast in (1.0, 1e5), height in (0.02, 2.0)
        state=(k2=(1e-8+1e-6im, contrast*im), mu=(1e-6,3e-6),
            sh=(1e-6im,0.1+1e-5im), k=E.outgoing_root.((1e-8+1e-6im,contrast*im)))
        angle=0.1
        weight=(height,separation=height/2,q=state.k[2])
        features=E.SpectralFeatures([1e-4,1.0,1e5])
        limit=max(10/height,4maximum(abs,state.k))
        for kind in (:Z,:phi,:voltage)
            kernel=E.EarthRadialSpectrum{kind,typeof(state),Float64}(state,0.3)
            original=E.earth_spectrum(Val(kind),Val(2),Val(2),state,height/2,height/2,0.3)
            transformed=E.RadializedEarthSpectrum(original,state.k[2],height)
            for k in (kernel,transformed)
                integral=E.SpectralIntegral(Val(:radial),k,weight,1.0;angle,features)
                bound=E.spectral_tail(integral)(limit)
                value,error=E.quadgk(t->abs(integral(cis(angle)*(limit+t/(1-t))))/(1-t)^2,
                    0.0,1.0;rtol=1e-7,atol=max(bound*1e-8,floatmin(Float64)))
                @test isfinite(bound) && value<=bound+error
            end
        end
    end
end

@testitem "Engine / analytic image tails include cosine Bessel and radial continuation" tags=[:unit] begin
    const E=LineCableModels.Engine
    angle=0.15
    amplitudes=ComplexF64[1+0.2im,-0.8+0.3im]
    poles=ComplexF64[0.3+0.4im,0.02-0.1im]
    for kind in (:cosine,:besselcosine,:radial)
        weight=kind===:radial ? (height=1.0,separation=0.5,q=0.7+0.2im) :
               kind===:besselcosine ? (height=1.0,separation=0.5,radius=0.1) :
               (height=1.0,separation=0.5)
        integral=E.SpectralIntegral(Val(kind),x->complex(exp(-x)),weight,1.0;
            angle,features=E.SpectralFeatures([1.0]))
        rotation=kind===:radial ? 1.0+0im : cis(angle)
        limit=8.0
        bound=E.cim_image_tail(integral,amplitudes,poles,rotation,1.0,limit)
        f=t->begin
            λ=cis(angle)*(limit+t/(1-t))
            abs(sum(E.cim_weighted_image(integral,a,b,λ,rotation) for (a,b) in zip(amplitudes,poles)))/(1-t)^2
        end
        value,error=E.quadgk(f,0.0,1.0;rtol=1e-8)
        @test value<=bound+error
        @test E.cim_image_tail(integral,amplitudes,poles,rotation,1.0,2limit)<bound
        @test isinf(E.cim_image_tail(integral,[1.0+0im],[-2.0+0im],rotation,1.0,limit))
    end
end

@testitem "Engine / finite coverage preserves inference and uncertain tails" tags=[:unit] begin
    using Measurements
    const E=LineCableModels.Engine
    state=(k2=(1e-8im,1.0im),k=E.outgoing_root.((1e-8im,1.0im)),
        mu=(1e-6,1e-6),sh=(1e-6im,0.1+1e-6im))
    function integral(state)
        kernel=E.earth_spectrum(Val(:Z),Val(2),Val(2),state,1.0,1.0,0.0)
        E.SpectralIntegral(Val(:cosine),kernel,(height=2.0,separation=1.0),1.0;
            features=E.SpectralFeatures([1e-4,1.0,1e5]))
    end
    bounded=@inferred integral(state)
    q=E.computation_options(E.SpectralIntegral,(method=:quad,options=(rtol=1e-10,))).options
    reference=E.integrate(Val(:quad),bounded,q,nothing)
    for tolerance in (1e-5,1e-8)
        controls=E.computation_options(E.SpectralIntegral,(method=:trapz,options=(rtol=tolerance,))).options
        estimate=@inferred E.spectral_estimate(Val(:trapz),bounded,controls,nothing)
        @test estimate.cutoff<1e5
        @test estimate.tail<=estimate.error<=tolerance*abs(estimate.value)
        @test abs(estimate.value-reference)<=estimate.error+1e-10*abs(reference)
    end
    uncertain=integral(merge(state,(k2=(state.k2[1],state.k2[2]*measurement(1.0,0.01)),)))
    @test E.spectral_tail(uncertain)===nothing
    @test @inferred(E.spectral_kernel_tail(uncertain.kernel,uncertain.weight,uncertain.angle))===nothing
end

@testitem "Engine / spectral manual construction budgets and adaptive defaults" tags=[:unit] begin
    const E=LineCableModels.Engine
    integral=E.SpectralIntegral(Val(:cosine),x->complex(exp(-x)),
        (height=1.0,separation=1.0),1.0)
    for method in (:trapz,:cim)
        automatic=E.computation_options(E.SpectralIntegral,(;method))
        @test automatic.options.samples===nothing
        manual=E.computation_options(E.SpectralIntegral,(;method,options=(samples=4096,)))
        estimate=E.spectral_estimate(manual.method,integral,manual.options,nothing)
        @test estimate.value≈0.4 rtol=1e-6
        @test 0<estimate.samples<=4096
        @test estimate.evaluations>=estimate.samples
        limited=E.computation_options(E.SpectralIntegral,(;method,options=(samples=16,)))
        @test_throws ErrorException E.spectral_estimate(limited.method,integral,limited.options,nothing)
    end
end

@testitem "Engine / projected complex pencils reuse rank without conjugating poles" tags=[:unit] begin
    const E=LineCableModels.Engine
    z=exp(-0.05+0.17im)
    data=z.^(0:64)
    pencil=E.cim_pencil(data,0.1,6.4,1e-9,16,nothing)
    integral=E.SpectralIntegral(Val(:cosine),x->complex(exp(-x)),
        (height=1.0,separation=0.0),1.0)
    for rank in (1,4,8)
        poles=E.cim_poles!(ComplexF64[],[pencil],rank,integral,1.0+0im)
        @test length(poles)==1
        @test only(poles)≈-log(z)/0.1 rtol=1e-12
    end
end
