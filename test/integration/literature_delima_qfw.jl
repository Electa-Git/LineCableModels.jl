@testitem "Engine / literature assimilation / De Lima image and prescribed qFW" begin
    using LineCableModels, SpecialFunctions, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance
    mu=4π*1e-7; epsilon=8.8541878128e-12
    for order in 0:3,z in (1e-10im,.5im,50im,-.5im,.01+1im,1.0+2im,10.0+0im)
        value=E._rotated_besselk(order,Complex{BigFloat}(z))
        @test value ≈ besselk(order,z) rtol=1e-12
        @test E.special_besselk(order,Complex{BigFloat}(z)) ≈ besselk(order,z) rtol=1e-12
    end
    function image_reference(s,h,r,kappa1,kappa2)
        g1=s*mu*kappa1; g2=s*mu*kappa2
        beta=sqrt(g2-g1); N=g2/g1; D=hypot(2h,r)
        S1=log(1+2/(beta*D))
        S2=2/(N+1)*log(1+(N+1)/(beta*D))
        S4=2log(2)+2N/(N+1)*log((1+(N+1)/(beta*D))/(1+2(N+1)/(beta*D)))
        base=log(2h/r)
        return (impedance=s*mu/(2π)*(base+S1-S2-S4),
            potential=s/(2π*kappa1)*(base-S4))
    end
    function spectral_reference(s,h,r,kx,kappa1,kappa2)
        g1=s*mu*kappa1; g2=s*mu*kappa2; N=g2/g1
        eta=sqrt(g1+kx^2); D=hypot(2h,r)
        Lambda=iszero(eta) ? log(D/r) : besselk(0,eta*r)-besselk(0,eta*D)
        S1,S2,S4=quadgk(0.,Inf;rtol=1e-11) do t
            lambda=t/h
            q1=sqrt(lambda^2+g1+kx^2); q2=sqrt(lambda^2+g2+kx^2)
            decay=exp(-h*q1); denominator=N*q1+q2
            weight=2cos(lambda*r)/h
            [decay^2/(q1+q2),decay^2/denominator,
                q2/q1*(decay-decay^2)/denominator]*weight
        end[1]
        return (impedance=s*mu/(2π)*(Lambda+S1+kx^2/g1*(S2+S4)),
            potential=s/(2π*kappa1)*(Lambda-S4))
    end
    for f in (50.,1e5,1e7),h in (.5,10.),rho in (10.,1000.),layer in (1,2)
        s=2π*f*im; r=.005
        args=([Inf,rho],[epsilon,10epsilon],[mu,mu],s)
        kappas=(s*epsilon,inv(rho)+s*10epsilon)
        k1,k2=layer==1 ? kappas : reverse(kappas)
        height=layer==1 ? h : -h
        pair=E.EarthPair(1,1,(height,height),r,(layer,layer))
        image_functor=EI.Formula(:DeLima2018)(args...,nothing)
        parameters=EI.earth_impedance(Val(:DeLima2018),Val(:parameters),image_functor,pair)
        reference=image_reference(s,h,r,k1,k2)
        @test parameters.impedance ≈ reference.impedance rtol=1e-10
        @test parameters.potential ≈ reference.potential rtol=1e-10
        @test image_functor(Val(:self),pair)==parameters.impedance
        @test EA.Formula(:DeLima2018)(args...,nothing)(Val(:self),pair)==parameters.potential
        negative=EI.Formula(:DeLima2018)(args[1:3]...,-s,nothing)
        @test negative(Val(:self),pair) ≈ conj(parameters.impedance) rtol=1e-12
        @test EA.Formula(:DeLima2018)(args[1:3]...,-s,nothing)(Val(:self),pair) ≈
            conj(parameters.potential) rtol=1e-12
        internal=E.InternalImpedance.Formula(:Schelkunoff1934)(0.,r,1.7241e-8,1.,s)(Val(:outer))
        wave=EI.propagation_constant(Val(:DeLima2018),Val(:image),parameters,internal,s)
        @test -wave.squared ≈ (internal+parameters.impedance)*s/parameters.potential rtol=1e-12
        @test real(-im*wave.Γ)>=0
        for kx in (wave.Γ,zero(s))
            reference=spectral_reference(s,h,r,kx,k1,k2)
            impedance=EI.Formula(:DeLima2018;approximation=:quasi_full_wave)
            potential=EA.Formula(:DeLima2018;approximation=:quasi_full_wave)
            actual=impedance(args...,kx)(Val(:self),pair)
            @test actual ≈ reference.impedance rtol=2e-7
            @test potential(args...,kx)(Val(:self),pair) ≈ reference.potential rtol=2e-7
            @test impedance(args[1:3]...,-s,-conj(kx))(Val(:self),pair) ≈ conj(actual) rtol=1e-9
            @test potential(args[1:3]...,-s,-conj(kx))(Val(:self),pair) ≈ conj(reference.potential) rtol=2e-7
        end
    end
    for T in (Float32,Float64,BigFloat),layer in (1,2),approximation in (:image,:quasi_full_wave)
        m=T(4)*T(π)/T(10)^7; e=T(epsilon); s=complex(zero(T),T(100)*T(π))
        args=(T[Inf,100],T[e,10e],T[m,m],s)
        height=T(layer==1 ? 10 : -10); r=T(.01)
        pair=E.EarthPair(1,1,(height,height),r,(layer,layer))
        for family in (EI,EA)
            value=family.Formula(:DeLima2018;approximation)(args...,zero(s))(Val(:self),pair)
            @test value isa Complex{T}
            @test isfinite(value)
        end
    end
    args=([Inf,100.],[epsilon,10epsilon],[mu,mu],100π*im)
    self=E.EarthPair(1,1,(1.,1.),.01,(1,1))
    mutual=E.EarthPair(1,2,(1.,1.),1.,(1,1))
    for family in (EI,EA)
        @test_throws ArgumentError family.Formula(:DeLima2018;approximation=:invalid)
        @test_throws ArgumentError family.Formula(:DeLima2018;approximation=:quasi_full_wave)(args...,nothing)
        @test_throws ArgumentError family.Formula(:DeLima2018)(args...,1.0+0im)
        @test_throws ArgumentError family.Formula(:DeLima2018)(args...,nothing)(Val(:mutual),mutual)
        @test_throws ArgumentError family.Formula(:DeLima2018)(args...,nothing)(Val(:self),mutual)
        @test_throws DomainError family.Formula(:DeLima2018)(args[1:2]...,[mu,2mu],args[4],nothing)
        @test_throws DomainError family.Formula(:DeLima2018)(args[1:3]...,zero(args[4]),nothing)
    end
    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    wire=build(CableDesign,"De Lima wire",Stack(Group(:core,Region(:core,Disk(.005),copper))))
    for height in (10.,-1.)
        system=build(LineCableSystem,[wire],[(0.,height)];connections=[Dict("core"=>1)])
        frequencies=[50.,1e5]; image_values=NamedTuple[]; waves=ComplexF64[]
        for f in frequencies
            s=2π*f*im
            pair=E.EarthPair(1,1,(height,height),.005,(height>0 ? 1 : 2,height>0 ? 1 : 2))
            image_functor=EI.Formula(:DeLima2018)(args[1:3]...,s,nothing)
            parameters=EI.earth_impedance(Val(:DeLima2018),Val(:parameters),image_functor,pair)
            push!(image_values,parameters)
            internal=E.InternalImpedance.Formula(:Schelkunoff1934)(0.,.005,1.7241e-8,1.,s)(Val(:outer))
            push!(waves,EI.propagation_constant(Val(:DeLima2018),Val(:image),parameters,internal,s).Γ)
        end
        for approximation in (:image,:quasi_full_wave)
            problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
                frequencies,temperature=20.0,Γ=approximation===:image ? nothing : waves)
            formulation=Formulation(internal_impedance=:Schelkunoff1934,
                earth_impedance=formula(:DeLima2018;approximation),
                earth_admittance=formula(:DeLima2018;approximation),
                options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
            result=compute(problem,formulation;options=(trace=true,))
            trace=details(result).trace
            for (k,f) in enumerate(frequencies)
                s=2π*f*im; kappas=(s*epsilon,.01+s*10epsilon)
                k1,k2=height>0 ? kappas : reverse(kappas)
                reference=approximation===:image ? image_values[k] :
                    spectral_reference(s,abs(height),.005,waves[k],k1,k2)
                internal=E.InternalImpedance.Formula(:Schelkunoff1934)(0.,.005,1.7241e-8,1.,s)(Val(:outer))
                @test trace.Zg[1,1,k] ≈ reference.impedance rtol=2e-7
                @test trace.Pg[1,1,k] ≈ reference.potential rtol=2e-7
                @test result.Z[1,1,k] ≈ internal+reference.impedance rtol=2e-7
                @test result.Y[1,1,k] ≈ s/reference.potential rtol=2e-7
            end
        end
    end
end
