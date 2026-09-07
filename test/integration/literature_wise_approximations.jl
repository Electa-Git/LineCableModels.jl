@testitem "Engine / literature assimilation / Wise analytical potential approximations" begin
    using LineCableModels, QuadGK, SpecialFunctions, LinearAlgebra
    E=LineCableModels.Engine; EA=E.EarthAdmittance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    function source_variables(f,H,x,er,sigma)
        omega=2π*f
        alpha=omega*μ0*sigma
        auxiliary=sqrt(1+im*omega*ε0*(er-1)/sigma)
        xi=abs(auxiliary); eta=angle(auxiliary)
        u=cis(eta+π/4)
        a=er-im*sigma/(omega*ε0)
        scale=xi*sqrt(alpha)
        return (;a,u,g=H*scale,left=(H-im*x)*scale,right=(H+im*x)*scale,scale)
    end
    function literal_nine(v)
        a,u=v.a,v.u
        r1=u/2*(1-sqrt(1-4/(a+1))); r2=u/2*(1+sqrt(1-4/(a+1)))
        li_scaled(z)=-exp(z)*expint(z)
        # Equation (9), followed by C=2(M+jN) from source (7).
        return sum(-r2*li_scaled(g*r1)+r1*li_scaled(g*r2)
            for g in (v.left,v.right))/((a+1)*(r2-r1))
    end
    function integral_reference(v,method,x)
        a,u,g=v.a,v.u,v.g
        transition=min(1.0,abs(u/a)*g)
        if method===:small_g
            return 2*(quadgk(t->inv(u+a*t),0.0,min(1.,abs(u/a)),1.0;rtol=1e-11)[1]+
                quadgk(t->exp(-g*t)/((1+a)*t),1.0,Inf;rtol=1e-11)[1])
        end
        return 2quadgk(0.0,transition,1.0,Inf;rtol=1e-11) do t
            lambda=t/g
            denominator=method===:integral ? sqrt(lambda^2+u^2)+a*lambda :
                method===:coarse ? u+a*lambda :
                ((a+1)*lambda^2+(a+1)*u*lambda+u*u)/(lambda+u)
            exp(-t)*cos(x*v.scale*lambda)/denominator/g
        end[1]
    end
    for (f,H,x,er,sigma) in ((1e3,10.,0.,15.,.01),(1e4,10.,3.,10.,.01),
            (1e4,10.,15.,10.,.01),(1e4,10.,100.,10.,.01),
            (1e6,1.,0.,5.,.001),(1e6,1.,10.,5.,.001),
            (1e4,1.,0.,3.,1e-8))
        v=source_variables(f,H,x,er,sigma)
        s=2π*f*im; args=([Inf,inv(sigma)],[ε0,er*ε0],[μ0,μ0],s,nothing)
        self=iszero(x); radius=.005
        heights=self ? (H/2,H/2) : (H/3,2H/3)
        pair=E.EarthPair(1,self ? 1 : 2,heights,self ? radius : x,(1,1))
        direct=self ? radius : hypot(x,heights[1]-heights[2])
        ideal=log(hypot(H,x)/direct)
        for method in (:integral,:rational,:coarse,:small_g)
            method in (:coarse,:small_g) && (!self || (method===:small_g && v.g>=1)) && continue
            recipe=EA.Formula(:Wise1948;approximation=method)
            value=recipe(args...)(self ? Val(:self) : Val(:mutual),pair)
            correction=integral_reference(v,method,x)
            @test value ≈ (ideal+correction)/(2π*ε0) rtol=(method===:integral ? 2e-8 : 1e-9)
            @test recipe(args[1:3]...,-s,nothing)(self ? Val(:self) : Val(:mutual),pair) ≈
                conj(value) rtol=1e-9
            if method===:rational && x<H/2
                @test correction ≈ literal_nine(v) rtol=1e-8
            end
            # Ordinary selected examples reproduce the paper's small correction,
            # without imposing its pointwise accuracy claim globally.
            @test isfinite(value)
        end
    end
    for T in (Float32,Float64,BigFloat),method in (:integral,:rational,:coarse,:small_g)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),T(2000)*T(π))
        args=(T[Inf,100],T[epsilon,15epsilon],T[mu,mu],s,nothing)
        pair=E.EarthPair(1,1,(T(5),T(5)),T(.005),(1,1))
        recipe=EA.Formula(:Wise1948;approximation=method)
        before=precision(BigFloat)
        value=recipe(args...)(Val(:self),pair)
        expected=(log(10/.005)+integral_reference(source_variables(1000.,10.,0.,15.,.01),
            method,0.))/(2π*ε0)
        @test value isa Complex{T}
        @test value ≈ expected rtol=(T===Float32 ? 2e-5 : 1e-9)
        @test precision(BigFloat)==before
    end
    # Continue the scaled E1 product across its apparent negative-axis cut.
    q=-1.0+1.0im
    for x in (.99,1.,1.01,3.),sense in (-1.,1.)
        z=complex(1.,sense*x); root=sense<0 ? conj(q) : q
        actual=EA.earth_potential_coefficient(Val(:Wise1948),Val(:laplace),root,z)
        expected=quadgk(t->exp(-z*t)/(t+root),0.,Inf;rtol=1e-12)[1]
        @test actual ≈ expected rtol=1e-10
    end
    # Nearly repeated quadratic roots use the same rational kernel.
    a=3.0+0im; beta=1.0+1im; H=2.; x=.3
    actual=EA.earth_potential_coefficient(Val(:Wise1948),Val(:analytical),a,beta,H,x,Val(:rational))
    expected=2quadgk(t->exp(-H*t)*cos(x*t)*(t+beta)/
        ((a+1)*t^2+(a+1)*beta*t+beta^2),0.,Inf;rtol=1e-12)[1]
    @test actual ≈ expected rtol=1e-10

    args=([Inf,100.],[ε0,15ε0],[μ0,μ0],2000π*im,nothing)
    mutual=E.EarthPair(1,2,(5.,5.),1.,(1,1))
    for method in (:coarse,:small_g)
        @test_throws ArgumentError EA.Formula(:Wise1948;approximation=method)(args...)(Val(:mutual),mutual)
    end
    @test_throws DomainError EA.Formula(:Wise1948;approximation=:small_g)(args...)(
        Val(:self),E.EarthPair(1,1,(10000.,10000.),.01,(1,1)))
    @test_throws ArgumentError EA.Formula(:Wise1948;approximation=:unknown)
    @test_throws DomainError EA.Formula(:Wise1948;approximation=:rational)(
        args[1:2]...,[μ0,2μ0],args[4:5]...)

    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    design=build(CableDesign,"Wise-wire",Group(:core,Region(:core,Disk(.005),copper)))
    for method in (:integral,:rational,:coarse,:small_g)
        # Vertically aligned distinct wires admit both zero-horizontal forms.
        system=build(LineCableSystem,[design,design],[(0.,5.),(0.,10.)];
            connections=[Dict("core"=>1),Dict("core"=>2)])
        problem=LineParametersProblem(system;earth_props=EarthModel(100.,15.,1.),
            frequencies=[50.,1000.],temperature=20.0)
        formulation=Formulation(earth_impedance=:Wise1934,
            earth_admittance=formula(:Wise1948;approximation=method),
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        result=compute(problem,formulation;options=(trace=true,))
        trace=details(result).trace
        @test size(result.Y)==(2,2,2)
        for k in 1:2
            s=2π*problem.frequencies[k]*im
            reference=[begin
                h1=(i==1 ? 5. : 10.); h2=(j==1 ? 5. : 10.)
                d=i==j ? .005 : abs(h1-h2)
                correction=integral_reference(source_variables(problem.frequencies[k],h1+h2,
                    0.,15.,.01),method,0.)
                (log((h1+h2)/d)+correction)/(2π*ε0)
            end for i in 1:2,j in 1:2]
            @test trace.Pg[:,:,k] ≈ reference rtol=1e-9
            @test result.Y[:,:,k] ≈ s*inv(reference) rtol=1e-9
            @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k]) rtol=1e-12
            @test minimum(eigvals(Symmetric(real.(result.Y[:,:,k]))))>=-1e-12
        end
    end
end
