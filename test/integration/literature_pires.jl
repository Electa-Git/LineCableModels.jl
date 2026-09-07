@testitem "Engine / literature assimilation / Pires double-exponential quadrature" begin
    using LineCableModels, SpecialFunctions, QuadGK, LinearAlgebra
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance
    for T in (Float32,Float64,BigFloat), a in (T(.01),one(T),T(100))
        tol=T===Float32 ? T(2e-5) : T(1e-10)
        for (f,expected) in ((x->exp(-a*x),inv(a)),
                (x->inv(one(T)+x)^2,one(T)),
                (x->exp(-x)/sqrt(x),sqrt(T(π))),
                (x->exp(-complex(a,a)*x),inv(complex(a,a))))
            result=E.double_exponential(f,T;rtol=tol)
            @test result.value ≈ expected rtol=10tol
            @test result.error>=0
            @test 3<=result.level<=12
        end
        for k in -3:3
            h=T(.25); t=T(k)*h
            x=tanh(T(π)/2*sinh(t))
            w=T(π)/2*h*cosh(t)/cosh(T(π)/2*sinh(t))^2
            lambda=(1+x)/(1-x)
            weight=2w/(1-x)^2
            @test lambda ≈ exp(T(π)*sinh(t)) rtol=10eps(T)
            @test weight ≈ h*T(π)*cosh(t)*lambda rtol=10eps(T)
        end
    end
    @test_throws ArgumentError E.double_exponential(exp,Float64;rtol=0.)
    @test_throws ArgumentError E.double_exponential(exp,Float64;maxlevel=2)
    @test_throws DomainError E.double_exponential(x->NaN,Float64)
    @test_throws ErrorException E.double_exponential(x->inv(1+x),Float64;maxlevel=6)

    μ0=4π*1e-7; ε0=8.8541878128e-12
    function reference(s,sigma,epsilon,hi,hj,x)
        kappa=sigma .+ s.*epsilon; g=s*μ0.*kappa
        H=hi+hj; d=hypot(x,hi-hj); D=hypot(x,H)
        function kernels(lambda,kind)
            q0,q1=sqrt.(lambda^2 .+ g)
            base=exp(-H*q1)*cos(lambda*x)
            return kind===:Z ? base/(q0+q1) :
                base*q0/(q1*(q0+g[1]/g[2]*q1))
        end
        # Source (5) and (6) independently evaluated, including lossy sea above seabed.
        direct=besselk(0,sqrt(g[2])*d)-besselk(0,sqrt(g[2])*D)
        spectral=map((:Z,:P)) do kind
            quadgk(lambda->kernels(lambda,kind),0.,1/H,Inf;rtol=1e-10)[1]
        end
        return s*μ0/(2π)*(direct+2spectral[1]),s/(2π*kappa[2])*(direct+2spectral[2])
    end
    for f in (10.,1e3,1e5,1e7), (sigma,epsilon) in
            (([0.,.01],[ε0,10ε0]),([5.,1.5],[81ε0,40ε0])),
            (hi,hj,x) in ((1.,1.,.07105),(1.,1.,1.),(.1,.4,.5),(10.,12.,2.))
        s=2π*f*im
        args=(inv.(sigma),epsilon,[μ0,μ0],s,nothing)
        self=hi==hj && x==.07105
        pair=E.EarthPair(1,self ? 1 : 2,(-hi,-hj),x,(2,2))
        kind=self ? Val(:self) : Val(:mutual)
        zref,pref=reference(s,sigma,epsilon,hi,hj,x)
        for (family,expected) in ((EI,zref),(EA,pref))
            adaptive=family.Formula(:Xue2018b)(args...)(kind,pair)
            de=family.Formula(:Xue2018b;quadrature=:double_exponential)(args...)(kind,pair)
            @test adaptive ≈ expected rtol=2e-7
            @test de ≈ expected rtol=2e-7
            @test de ≈ adaptive rtol=2e-7
            negative=family.Formula(:Xue2018b;quadrature=:double_exponential)(
                args[1:3]...,-s,nothing)(kind,pair)
            @test negative ≈ conj(de) rtol=2e-7
        end
    end
    for T in (Float32,Float64,BigFloat),family in (EI,EA)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),T(100)*T(π))
        args=(T[Inf,100],T[epsilon,10epsilon],T[mu,mu],s,nothing)
        pair=E.EarthPair(1,1,(-one(T),-one(T)),T(.07105),(2,2))
        de=family.Formula(:Xue2018b;quadrature=:double_exponential)(args...)(Val(:self),pair)
        adaptive=family.Formula(:Xue2018b)(args...)(Val(:self),pair)
        @test de isa Complex{T}
        @test de ≈ adaptive rtol=(T===Float32 ? 2e-4 : 2e-7)
        @test_throws ArgumentError family.Formula(:Xue2018b;quadrature=:unknown)
    end

    copper=Material(:conductor,1.7e-8,1.,1.,20.,0.)
    sheath=Material(:conductor,21e-8,1.,1.,20.,0.)
    inner=Material(:insulator,Inf,3.5,1.,20.,0.)
    outer=Material(:insulator,Inf,8.,1.,20.,0.)
    design=build(CableDesign,"Pires-cable",Stack(
        Group(:core,Region(:core,Disk(.03395),copper)),
        Region(:insulation,Shell(.06065-.03395),inner),
        Group(:sheath,Region(:sheath,Shell(.06465-.06065),sheath)),
        Region(:jacket,Shell(.07105-.06465),outer)))
    system=build(LineCableSystem,[design,design],[(0.,-1.),(1.,-1.)];
        connections=[Dict("core"=>1,"sheath"=>2),Dict("core"=>3,"sheath"=>4)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
        frequencies=[50.,1e5],temperature=20.)
    results=map((:adaptive,:double_exponential)) do quadrature
        formulation=Formulation(earth_impedance=formula(:Xue2018b;quadrature),
            earth_admittance=formula(:Xue2018b;quadrature),
            insulation_admittance=:Ametani1980,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        compute(problem,formulation;options=(trace=true,))
    end
    @test results[1].Z ≈ results[2].Z rtol=2e-7
    @test results[1].Y ≈ results[2].Y rtol=2e-7
    trace=details(results[2]).trace
    for k in 1:2
        s=2π*problem.frequencies[k]*im
        @test results[2].Y[:,:,k] ≈ s*inv(trace.P[:,:,k]) rtol=1e-10
        @test results[2].Z[:,:,k] ≈ transpose(results[2].Z[:,:,k])
        @test minimum(eigvals(Symmetric(real.(results[2].Z[:,:,k]))))>=-1e-12
        @test minimum(eigvals(Symmetric(real.(results[2].Y[:,:,k]))))>=-1e-12
    end
end
