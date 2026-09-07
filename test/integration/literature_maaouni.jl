@testitem "Engine / literature assimilation / Maaouni overhead images and potential" begin
    using LineCableModels, SpecialFunctions, QuadGK, LinearAlgebra
    E=LineCableModels.Engine; EA=E.EarthAdmittance; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    Q(w)=exp(-w)*expint(-w)
    function printed_G(N,k0,H,x)
        A=sqrt(1-N); b=-im/sqrt(1+N)
        z=k0*complex(H,x); zbar=conj(z)
        P(b,z)=-(1-N)/(2b)*(log(1+2/(z*A))+Q(b*z+2b/A))+
            inv(z)+b*Q(b*z)*(1+(1-N)/(2b*b))
        return N/(2*(N*N-1))*(Q(b*z)+Q(b*zbar))-
            (P(b,z)+P(b,zbar)-P(-b,z)-P(-b,zbar)-
             N*b*(Q(-b*z)+Q(-b*zbar)))/(2b*(N*N-1))
    end
    function approximated_integral(N,beta,H,x)
        return 2quadgk(0.0,Inf;rtol=1e-11) do t
            lambda=t/H
            replacement=iszero(lambda) ? beta :
                lambda+beta^2*(-expm1(-2lambda/beta))/(2lambda)
            numerator=N*lambda-replacement
            denominator=(N*N-1)*lambda^2-beta^2
            exp(-t)*cos(x*lambda)*numerator/denominator/H
        end[1]
    end
    for (frequency,x,H,er,sigma) in ((1e3,0.,10.,15.,.01),
            (1e4,5.,10.,10.,.01),(1e4,7.,10.,10.,.001),
            (1e6,8.,20.,5.,.01),(1e4,0.,.2,10.,.001),(1e6,0.,.1,10.,.01))
        s=2π*frequency*im
        args=([Inf,inv(sigma)],[ε0,er*ε0],[μ0,μ0],s,nothing)
        f=EA.Formula(:Maaouni2001)(args...)
        N=er+sigma/(s*ε0); k0=imag(s)*sqrt(μ0*ε0)
        beta=k0*sqrt(1-N)
        self=iszero(x)
        heights=self ? (H/2,H/2) : (H/3,2H/3)
        radius=min(.001,H/20)
        pair=E.EarthPair(1,self ? 1 : 2,heights,self ? radius : x,(1,1))
        kind=self ? Val(:self) : Val(:mutual)
        direct=self ? radius : hypot(x,heights[1]-heights[2])
        ideal=log(hypot(H,x)/direct)
        literal=printed_G(N,k0,H,x)
        integral=approximated_integral(N,beta,H,x)
        kernel=EA.earth_potential_coefficient(Val(:Maaouni2001),Val(:kernel),N,k0,H,x,1.0)
        @test kernel ≈ literal rtol=1e-8 atol=1e-14
        @test kernel ≈ integral rtol=2e-8 atol=1e-14
        @test f(kind,pair) ≈ (ideal+literal)/(2π*ε0) rtol=1e-12
        @test imag(f(kind,pair)) >= 0
        @test EA.Γ(f) ≈ k0 rtol=1e-12
        explicit=EA.Formula(:Maaouni2001)(args[1:4]...,complex(k0))
        @test explicit(kind,pair) ≈ f(kind,pair) rtol=1e-12
        negative=EA.Formula(:Maaouni2001)(args[1:3]...,-s,nothing)
        @test negative(kind,pair) ≈ conj(f(kind,pair)) rtol=1e-12
        # Source (7), without the approximation (23), is the Wise integral.
        parent=2quadgk(0.0,Inf;rtol=1e-11) do t
            lambda=t/H
            exp(-t)*cos(x*lambda)/(N*lambda+sqrt(lambda^2+beta^2))/H
        end[1]
        @test abs(kernel-parent)/abs(parent)<.005
        if !self
            @test EA.Formula(:Wise1948)(args...)(Val(:mutual),pair) ≈
                (ideal+parent)/(2π*ε0) rtol=1e-9
        end
        # Maaouni (12) and Ametani (17)–(18), using the same SI constants.
        image=sqrt(x*x+(H+2/beta)^2)
        J=log(image/hypot(H,x))
        expected=s*μ0/(2π)*(ideal+J)
        @test EI.Formula(:Pettersson1994)(args...)(kind,pair) ≈ expected rtol=1e-12
        @test EI.Formula(:Ametani2014)(args...)(kind,pair) ≈ expected rtol=1e-12
    end
    for T in (Float32,Float64,BigFloat)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),T(2000)*T(π))
        f=EA.Formula(:Maaouni2001)(T[Inf,100],T[epsilon,15epsilon],T[mu,mu],s,nothing)
        pair=E.EarthPair(1,1,(T(5),T(5)),T(1)/100,(1,1))
        before=precision(BigFloat)
        value=f(Val(:self),pair)
        @test value isa Complex{T}
        @test isfinite(value)
        @test precision(BigFloat)==before
        for w in (complex(T(1)/10,T(1)/5),complex(T(3),T(-2)),complex(T(-100),T(20)))
            reference=Q(ComplexF64(w))
            @test E.scaled_expint_negative(w) ≈ reference rtol=(T===Float32 ? 2e-5 : 3e-12)
        end
        w=complex(T(-600),T(300))
        asymptote=-sum(T(factorial(big(k)))/w^(k+1) for k in 0:12)
        @test E.scaled_expint_negative(w) ≈ asymptote rtol=(T===Float32 ? 2e-5 : 1e-12)
    end
    s=2π*1e6im
    args=([Inf,100.],[ε0,10ε0],[μ0,μ0],s,nothing)
    f=EA.Formula(:Maaouni2001)(args...)
    @test_throws DomainError f(Val(:mutual),E.EarthPair(1,2,(1.,1.),100.,(1,1)))
    @test_throws ArgumentError f(Val(:mutual),E.EarthPair(1,2,(1.,-1.),.1,(1,2)))
    @test_throws ArgumentError EA.Formula(:Maaouni2001)(args[1:4]...,0.0im)
    @test_throws DomainError EA.Formula(:Maaouni2001)([100.,100.],args[2:5]...)
    @test_throws DomainError EA.Formula(:Maaouni2001)(args[1:2]...,[μ0,2μ0],args[4:5]...)

    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    design=build(CableDesign,"Maaouni-wire",Group(:core,Region(:core,Disk(.01),copper)))
    system=build(LineCableSystem,[design,design],[(0.,10.),(3.,15.)];
        connections=[Dict("core"=>1),Dict("core"=>2)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
        frequencies=[50.,1e5],temperature=20.0)
    formulation=Formulation(earth_impedance=:Pettersson1994,earth_admittance=:Maaouni2001,
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    result=compute(problem,formulation;options=(trace=true,))
    trace=details(result).trace
    @test size(result.Y)==(2,2,2)
    for k in 1:2
        @test all(isfinite,result.Y[:,:,k])
        @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k]) rtol=1e-12
        @test result.Y[:,:,k] ≈ 2π*problem.frequencies[k]*im*inv(trace.P[:,:,k]) rtol=1e-12
        @test minimum(eigvals(Symmetric(real.(result.Y[:,:,k])))) >= -1e-12
    end
end
