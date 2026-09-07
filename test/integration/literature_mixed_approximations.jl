@testitem "Engine / literature assimilation / dimensional mixed approximations" begin
    using LineCableModels, SpecialFunctions, LinearAlgebra
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    for T in (Float32,Float64,BigFloat),frequency in (1,50,10000),
            (ha,hg,x) in ((10,1,0),(2,5,3),(1,1,7))
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),2T(π)*T(frequency)); g=sqrt(s*mu/T(100))
        h1,h2,spacing=T(ha),T(hg),T(x); d=hypot(spacing,h1+h2)
        args=(T[Inf,100],T[epsilon,10epsilon],T[mu,mu],s,nothing)
        pair=E.EarthPair(1,2,(h1,-h2),spacing,(1,2))
        reverse=E.EarthPair(2,1,(-h2,h1),spacing,(2,1))
        tolerance=T===Float32 ? T(2e-6) : T(1e-12)
        for approximation in (:ccitt,:wedepohl)
            recipe=EI.Formula(:Uribe2008;approximation)
            value=recipe(args...)(Val(:mutual),pair)
            if approximation===:ccitt
                ke=complex(zero(T),-one(T))*g
                reference=s*mu/(2T(π))*(log(T(1851)/1000/(im*ke*d))+
                    2im*ke*(h1-h2)/3)
            else
                # The dimensional (6e) and the parent Bessel Euler factor.
                p=inv(g); euler_factor=exp(T(Base.MathConstants.eulergamma))
                reference=s*mu/(2T(π))*(-log(euler_factor*d/(2p))+T(.5)-
                    4*(h1+h2)/(3p))
            end
            @test value isa Complex{T}
            @test value ≈ reference rtol=tolerance
            @test recipe(args...)(Val(:mutual),reverse) ≈ reference rtol=tolerance
            @test recipe(args[1:3]...,-s,nothing)(Val(:mutual),pair) ≈ conj(reference) rtol=tolerance
            if T===Float64
                exact=EI.Formula(:Pollaczek1926)(args...)(Val(:mutual),pair)
                relative_error=abs(value/exact-1)
                @test isfinite(relative_error)
                frequency<=50 && @test relative_error<.04
            end
            # Exact same-medium leaves are retained, not reimplemented.
            for (heights,layers) in (((h1,h1),(1,1)),((-h2,-h2),(2,2)))
                self=E.EarthPair(1,1,heights,T(.01),layers)
                @test recipe(args...)(Val(:self),self) ≈
                    EI.Formula(:Pollaczek1926)(args...)(Val(:self),self) rtol=tolerance
            end
        end
        approximate=EI.Formula(:Ametani2009;approximation=:power_frequency)
        value=approximate(args...)(Val(:mutual),pair)
        reference=s*mu/(2T(π))*(log(2/(abs(g)*d))-complex(zero(T),T(π)/4))
        @test value isa Complex{T}
        @test value ≈ reference rtol=tolerance
        @test approximate(args...)(Val(:mutual),reverse) ≈ reference rtol=tolerance
        @test approximate(args[1:3]...,-s,nothing)(Val(:mutual),pair) ≈ conj(reference) rtol=tolerance
        parent=EI.Formula(:Ametani2009)(args...)(Val(:mutual),pair)
        image=sqrt((h1+h2+2/g)^2+spacing^2)
        @test parent ≈ s*mu/(2T(π))*exp(-h2*g)*log(image/d) rtol=tolerance
        @test EI.Formula(:Ametani2009)(args...)(Val(:mutual),reverse) ≈ parent rtol=tolerance
    end

    # The low-frequency limit is taken from the full exponential image,
    # and the dimensional constants agree with the rounded published form.
    pair=E.EarthPair(1,2,(10.,-1.),3.,(1,2)); d=hypot(11.,3.)
    errors=Float64[]
    for f in (1e-4,1e-6,1e-8)
        s=2π*f*im; args=([Inf,100.],[ε0,10ε0],[μ0,μ0],s,nothing)
        leading=EI.Formula(:Ametani2009;approximation=:power_frequency)(args...)(Val(:mutual),pair)
        parent=EI.Formula(:Ametani2009)(args...)(Val(:mutual),pair)
        push!(errors,abs(leading/parent-1))
        printed=f*(1+im*(8.253+.628*log(100/(f*d*d))))*1e-6
        @test real(leading) ≈ real(printed) rtol=.014
        @test imag(leading) ≈ imag(printed) rtol=.001
        # Independent Bessel expansion fixes the multiplicative Euler factor.
        g=sqrt(s*μ0/100)
        reference=s*μ0/(2π)*(besselk(0,g*d)+.5-4g*11/3)
        wed=EI.Formula(:Uribe2008;approximation=:wedepohl)(args...)(Val(:mutual),pair)
        @test wed ≈ reference rtol=1e-8
    end
    @test errors[2]<errors[1]/5
    @test errors[3]<errors[2]/5
    @test_throws ArgumentError EI.Formula(:Uribe2008;approximation=:unknown)
    @test_throws ArgumentError EI.Formula(:Ametani2009;approximation=:unknown)
    args=([Inf,100.],[ε0,10ε0],[μ0,μ0],100π*im,nothing)
    @test_throws ArgumentError EI.Formula(:Uribe2008)(args[1:4]...,1.0+0im)
    @test_throws DomainError EI.Formula(:Uribe2008)(args[1:2]...,[μ0,2μ0],args[4:5]...)

    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    insulation=Material(:insulator,Inf,2.3,1.0,20.0,0.0)
    design=build(CableDesign,"mixed-wire",Stack(
        Group(:core,Region(:core,Disk(.005),copper)),Region(:coat,Shell(.003),insulation)))
    system=build(LineCableSystem,[design,design],[(0.,10.),(3.,-1.)];
        connections=[Dict("core"=>1),Dict("core"=>2)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
        frequencies=[1.,50.],temperature=20.0)
    for (identifier,approximation) in ((:Uribe2008,:ccitt),(:Uribe2008,:wedepohl),
            (:Ametani2009,:image),(:Ametani2009,:power_frequency))
        formulation=Formulation(earth_impedance=formula(identifier;approximation),
            earth_admittance=:IdealGround,insulation_admittance=:Ametani1980,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        result=compute(problem,formulation;options=(trace=true,))
        trace=details(result).trace
        @test size(result.Z)==(2,2,2)
        for k in 1:2
            s=2π*problem.frequencies[k]*im
            recipe=EI.Formula(identifier;approximation)
            f=recipe([Inf,100.],[ε0,10ε0],[μ0,μ0],s,nothing)
            @test trace.Zg[1,2,k] ≈ f(Val(:mutual),pair) rtol=1e-10
            @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k]) rtol=1e-12
            @test result.Z[:,:,k]-trace.Zin[:,:,k] ≈ trace.Zg[:,:,k] rtol=1e-10
            @test minimum(eigvals(Symmetric(real.(result.Z[:,:,k]))))>=0
        end
    end
end
