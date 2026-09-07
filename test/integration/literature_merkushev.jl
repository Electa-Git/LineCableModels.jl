@testitem "Engine / literature assimilation / Merkushev linear ACSR" begin
    using LineCableModels,SpecialFunctions,QuadGK,LinearAlgebra
    E=LineCableModels.Engine;II=E.InternalImpedance
    # Integrate the same azimuthal conductivity over a unit strand disk
    # centred at radius 2, instead of the source's radial-sector integration.
    form_factor(theta)=quadgk(t->2t/sqrt(
        (1+theta^2*(4+t^2))^2-(4theta^2*t)^2),0.,1.;rtol=1e-12)[1]
    function source_reference(R,h,rhos,rhoa,mur,s)
        theta=2π*R/h
        Sz=6π*R^2/rhoa*form_factor(theta);Sc=π*R^2/rhos
        if iszero(s)
            return complex(inv(Sc+Sz))
        end
        k=sqrt(-s*4π*1e-7*mur/rhos)
        j0=besselj(0,k*R);j1=besselj(1,k*R)
        if isinf(h)
            steel=rhos*k/(2π*R)*j0/j1
            return inv(inv(steel)+Sz)
        end
        A=k*R*j0/j1
        B=2Sc/(theta*Sz)-theta*k*R*j1/j0
        gamma=B/A
        k*rhos/h*gamma/(1+gamma*theta)*j0/j1
    end
    for T in (Float32,Float64,BigFloat),f in (0.,.01,50.,250.),
            mur in (1.,100.,10000.),pitch in (.04,.09,Inf)
        R=T(.0019);h=T(pitch);rhos=T(1/7.3e6);rhoa=T(1/3.6e7)
        s=complex(zero(T),T(2)*T(π)*T(f))
        functor=II.Formula(:Merkushev2015)(R,h,rhos,rhoa,T(mur),s)
        value=functor(Val(:outer))
        reference=source_reference(Float64(R),Float64(h),Float64(rhos),
            Float64(rhoa),mur,ComplexF64(s))
        tol=T===Float32 ? 4e-5 : 4e-11
        @test value ≈ reference rtol=tol
        @test value isa Complex{T}
        @test real(value)>0 && imag(value)>=0
        @test functor.state.Q ≈ form_factor(2π*Float64(R)/Float64(h)) rtol=tol
        @test iszero(functor(Val(:inner))) && iszero(functor(Val(:mutual)))
        negative=II.Formula(:Merkushev2015)(R,h,rhos,rhoa,T(mur),-s)
        @test negative(Val(:outer)) ≈ conj(value) rtol=tol
    end
    for h in (.04,.09)
        z=II.Formula(:Merkushev2015)(.0019,h,1/7.3e6,1/3.6e7,100.,100π*im)(Val(:outer))
        scaled=II.Formula(:Merkushev2015)(.0038,2h,1/7.3e6,1/3.6e7,100.,25π*im)(Val(:outer))
        @test scaled ≈ z/4 rtol=1e-11
    end
    @test_throws DomainError II.Formula(:Merkushev2015)(
        .0019,.09,1/7.3e6,1/3.6e7,100.,2e6π*im)
    @test_throws ArgumentError II.Formula(:Merkushev2015)(0.,.0057,1e-8,1.,100π*im)

    function acsr_design(mur;radius=.0019,pitch=.09,wirecount=6,hand=1)
        steel=Material(:conductor,1/7.3e6,1.,mur,20.,.002)
        aluminium=Material(:conductor,1/3.6e7,1.,1.,20.,.004)
        path=isinf(pitch) ? nothing : Helix(Pitch(pitch);dir=hand)
        design=build(CableDesign,"linear-acsr",Group(:core,Stack(
            Region(:steel,Disk(radius),steel),
            Group(:aluminium,Region(:strand,Disk(radius),aluminium);
                pattern=Ring(wirecount;r=2radius),path=path))))
        return design
    end
    unsupported=acsr_design(100.;wirecount=5)
    @test only(E.flatten(LineCableModelsCoaxial(),unsupported).conductors).acsr===nothing
    for pitch in (.09,Inf),hand in (-1,1),temperature in (20.,80.),
            correction in (false,true)
        design=acsr_design(100.;pitch,hand)
        blueprint=E.flatten(LineCableModelsCoaxial(),design)
        @test length(blueprint.conductors)==1
        profile=only(blueprint.conductors).acsr
        @test profile!==nothing
        @test profile.radius≈.0019
        @test profile.pitch==pitch
        system=build(LineCableSystem,[design,design],[(-1.,10.),(1.,10.)];
            connections=[Dict("core"=>1),Dict("core"=>2)])
        problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
            frequencies=[1.,50.,250.],temperature)
        formulation=Formulation(internal_impedance=:Merkushev2015,
            earth_impedance=:Carson1926,earth_admittance=:Wise1948,
            options=(temperature_correction=correction,reduce_bundle=false,
                kron_reduction=false,ideal_transposition=false))
        result=compute(problem,formulation;options=(trace=true,))
        trace=details(result).trace
        for k in 1:3
            rho1=profile.core.rho*(correction ? 1+profile.core.alpha*(temperature-20) : 1)
            rho2=profile.strands.rho*(correction ? 1+profile.strands.alpha*(temperature-20) : 1)
            expected=II.Formula(:Merkushev2015)(.0019,pitch,rho1,rho2,100.,
                2π*problem.frequencies[k]*im)(Val(:outer))
            @test trace.Zin[:,:,k] ≈ expected*Matrix{Float64}(I,2,2) rtol=1e-11
            @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k])
            @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k])
        end
    end
end
