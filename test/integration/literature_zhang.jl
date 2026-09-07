@testitem "Engine / literature assimilation / Zhang buried-wire potential" begin
    using LineCableModels, SpecialFunctions, QuadGK, LinearAlgebra
    E=LineCableModels.Engine; EA=E.EarthAdmittance; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    function original_reference(s,h,r,sigma,epsilon;extracted=false)
        g0=s*s*μ0*ε0; g1=s*μ0*(sigma+s*epsilon)
        kx1=-s*s*μ0*epsilon; a=s*μ0*sigma
        g=sqrt(a); threshold=10abs(g); transition=abs(g0/g1*g)
        function kernels(lambda;subtract_direct=false)
            u0=lambda; u1=sqrt(lambda^2+a); decay=exp(-2h*u1)
            F=(subtract_direct ? -decay : -expm1(-2h*u1))/u1+
                2decay/(u1+u0)
            G=2u1*(g1-g0)*decay/((u1+u0)*(u1*g0+u0*g1))
            return (F+G)*cos(lambda*r)
        end
        if extracted
            finite=quadgk(kernels,0.0,min(transition,threshold),threshold;rtol=1e-11)[1]
            moment(n,z)=threshold^(1-n)*real(expint(n,threshold*z))
            tail=moment(1,r*im)-a/2*moment(3,r*im)+a/2*moment(3,2h+r*im)+
                (g1-g0)/(g1+g0)*(moment(1,2h+r*im)+(g1-kx1)/2*moment(3,2h+r*im))
            integral=finite+tail
        else
            correction=quadgk(lambda->kernels(lambda;subtract_direct=true),
                0.0,min(transition,1/(2h)),1/(2h),Inf;rtol=1e-11)[1]
            integral=besselk(0,g*r)+correction
        end
        return s/(2π*(sigma+s*epsilon))*integral
    end
    for f in (1.,50.,1e5,1e7),h in (.05,1.,20.),sigma in (.0001,.01),
            evaluation in (:integral,:asymptotic_tail)
        s=2π*f*im; r=.005; epsilon=10ε0
        args=([Inf,inv(sigma)],[ε0,epsilon],[μ0,μ0],s,nothing)
        pair=E.EarthPair(1,1,(-h,-h),r,(2,2))
        recipe=EA.Formula(:Zhang2017;evaluation)
        functor=recipe(args...)
        coefficient=EA.earth_potential_coefficient(Val(:Zhang2017),Val(:coefficient),functor,pair)
        reference=original_reference(s,h,r,sigma,epsilon;
            extracted=evaluation===:asymptotic_tail)
        @test coefficient ≈ reference rtol=1e-8
        @test isfinite(coefficient)
        negative=recipe(args[1:3]...,-s,nothing)
        @test EA.earth_potential_coefficient(Val(:Zhang2017),Val(:coefficient),
            negative,pair) ≈ conj(coefficient) rtol=1e-10
        if real(s/reference)<-1e-10
            @test_throws DomainError functor(Val(:self),pair)
            @test_throws DomainError negative(Val(:self),pair)
        else
            @test functor(Val(:self),pair) ≈ reference rtol=1e-8
            @test real(s/functor(Val(:self),pair))>=-1e-10
        end
        # Source (17) and (18) are existing logarithm and scalar conversion.
        z=EI.Formula(:Petrache2005)(args...)(Val(:self),pair)
        zref=s*μ0/(2π)*log((1+sqrt(s*μ0*(sigma+s*epsilon))*r)/
            (sqrt(s*μ0*(sigma+s*epsilon))*r))
        @test z ≈ zref rtol=1e-12
        support=(functor,p)->EI.earth_impedance(Val(:Petrache2005),Val(:self),functor,p)
        p=EA.Formula(:Vance1978;impedance=support)(args...)(Val(:self),pair)
        @test s/p ≈ s*μ0*(sigma+s*epsilon)/zref rtol=1e-12
    end
    for T in (Float32,Float64,BigFloat),evaluation in (:integral,:asymptotic_tail)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),T(100)*T(π))
        args=(T[Inf,100],T[epsilon,10epsilon],T[mu,mu],s,nothing)
        pair=E.EarthPair(1,1,(-one(T),-one(T)),T(.01),(2,2))
        recipe=EA.Formula(:Zhang2017;evaluation)
        value=recipe(args...)(Val(:self),pair)
        @test value isa Complex{T}
        @test value ≈ original_reference(ComplexF64(s),1.,.01,.01,Float64(10epsilon);
            extracted=evaluation===:asymptotic_tail) rtol=(T===Float32 ? 5e-5 : 1e-8)
        for n in (1,3),z in (complex(T(.1),T(.2)),complex(T(3),T(.2)),
                complex(T(0),T(30)),complex(T(100),T(50)))
            moment=EA.earth_potential_coefficient(Val(:Zhang2017),Val(:tail_moment),n,z)
            @test moment ≈ expint(n,ComplexF64(z)) rtol=(T===Float32 ? 2e-5 : 1e-10) atol=1e-40
        end
    end
    args=([Inf,100.],[ε0,10ε0],[μ0,μ0],100π*im,nothing)
    recipe=EA.Formula(:Zhang2017)
    for pair in (E.EarthPair(1,2,(-1.,-1.),1.,(2,2)),
            E.EarthPair(1,2,(1.,-1.),1.,(1,2)))
        @test_throws ArgumentError recipe(args...)(Val(:self),pair)
        @test_throws ArgumentError recipe(args...)(Val(:mutual),pair)
    end
    @test_throws ArgumentError EA.Formula(:Zhang2017;evaluation=:unknown)
    @test_throws DomainError EA.Formula(:Zhang2017;source_radius=0)
    @test_throws ArgumentError recipe(args[1:4]...,1.0+0im)
    @test_throws DomainError recipe(args[1:2]...,[μ0,2μ0],args[4:5]...)
    @test_throws DomainError recipe([Inf,Inf],args[2:5]...)

    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    insulation=Material(:insulator,Inf,2.3,1.0,20.0,0.0)
    design=build(CableDesign,"Zhang-wire",Stack(
        Group(:core,Region(:core,Disk(.005),copper)),Region(:coat,Shell(.003),insulation)))
    system=build(LineCableSystem,[design],[(0.,-1.)];connections=[Dict("core"=>1)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
        frequencies=[50.,1e5],temperature=20.0)
    for evaluation in (:integral,:asymptotic_tail),radius in (nothing,.005)
        formulation=Formulation(earth_impedance=:Petrache2005,
            earth_admittance=formula(:Zhang2017;evaluation,source_radius=radius),
            insulation_admittance=:Ametani1980,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        result=compute(problem,formulation;options=(trace=true,))
        trace=details(result).trace
        @test size(result.Y)==(1,1,2)
        for k in 1:2
            s=2π*problem.frequencies[k]*im
            source_r=radius===nothing ? .008 : radius
            expected=original_reference(s,1.,source_r,.01,10ε0;
                extracted=evaluation===:asymptotic_tail)
            @test trace.Pg[1,1,k] ≈ expected rtol=1e-8
            @test trace.Pin[1,1,k] ≈ log(.008/.005)/(2π*2.3ε0) rtol=1e-12
            @test result.Y[1,1,k] ≈ s/(expected+trace.Pin[1,1,k]) rtol=1e-8
        end
    end
end
