@testitem "Engine / literature assimilation / Papadopoulos homogeneous F and G" begin
    using LineCableModels,QuadGK,SpecialFunctions,LinearAlgebra
    E=LineCableModels.Engine;EI=E.EarthImpedance;EA=E.EarthAdmittance
    μ0=4π*1e-7;ε0=8.8541878128e-12
    function source_integrals(s,mu,eps,rho,kx,pair)
        g=s.*mu.*(inv.(rho).+s.*eps)
        # Use the frequency-prescribed side of the propagating branch.
        root=z->iszero(imag(z))&&real(z)<0 ?
            complex(0.,sign(imag(s))*sqrt(-real(z))) : sqrt(z)
        h1,h2=abs.(pair.heights);H=h1+h2;x=pair.separation
        a1zero=root(g[2]+kx^2)
        base=besselk(0,a1zero*hypot(x,h1-h2))-
            besselk(0,a1zero*hypot(x,H))
        function fields(t)
            a0=root(t^2+g[1]+kx^2);a1=root(t^2+g[2]+kx^2)
            A=a1*mu[1]+a0*mu[2]
            B=a1*g[1]*mu[2]+a0*g[2]*mu[1]
            image=exp(-a1*H)*cos(x*t)
            F=2mu[1]/A*image
            G=2mu[1]*mu[2]*a1*(g[2]-g[1])/(A*B)*image
            [F,G]
        end
        point=iszero(imag(g[1]+kx^2))&&real(g[1]+kx^2)<0 ?
            sqrt(-real(g[1]+kx^2)) : 0.
        bounds=point>0 ? (0.,point,Inf) : (0.,Inf)
        integrals=quadgk(fields,bounds...;rtol=1e-10)[1]
        return (Z=s*mu[2]/(2π)*(base+integrals[1]),
            P=s/(2π*(1/rho[2]+s*eps[2]))*(base+sum(integrals)))
    end
    for frequency in (50.,1e5,1e7),mur in (.3,1.,5.),multiple in (0.,1.,2.)
        s=2π*frequency*im
        rho=[Inf,100.];epsilon=[ε0,10ε0];mu=[μ0,mur*μ0]
        kx=multiple*abs(s)*sqrt(mu[2]*epsilon[2])
        gamma=complex(kx)
        for pair in (E.EarthPair(1,2,(-1.,-3.),.7,(2,2)),
                E.EarthPair(1,1,(-1.,-1.),.02,(2,2)))
            expected=source_integrals(s,mu,epsilon,rho,kx,pair)
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            for (family,quantity) in ((EI,:Z),(EA,:P))
                leaf=family.Formula(:Papadopoulos2010b)(rho,epsilon,mu,s,gamma)
                value=leaf(kind,pair)
                @test value ≈ getproperty(expected,quantity) rtol=8e-7
                negative=family.Formula(:Papadopoulos2010b)(rho,epsilon,mu,-s,gamma)
                @test negative(kind,pair) ≈ conj(value) rtol=8e-7
                if multiple==1
                    default=family.Formula(:Papadopoulos2010b)(rho,epsilon,mu,s,nothing)
                    @test default(kind,pair) ≈ value rtol=8e-7
                end
            end
        end
    end
    for T in (Float32,Float64,BigFloat),family in (EI,EA)
        args=(T[Inf,100],T[ε0,10ε0],T[μ0,3μ0],
            complex(zero(T),T(100)*T(π)),nothing)
        pair=E.EarthPair(1,2,(T(-1),T(-3)),T(.7),(2,2))
        value=family.Formula(:Papadopoulos2010b)(args...)(Val(:mutual),pair)
        @test value isa Complex{T}
        @test isfinite(value)
    end
    copper=Material(:conductor,1.7241e-8,1.,1.,20.,0.)
    insulation=Material(:insulator,Inf,2.3,1.,20.,0.)
    design=build(CableDesign,"permeable-ground",Stack(
        Group(:core,Region(:core,Disk(.01),copper)),
        Region(:coat,Shell(.002),insulation)))
    system=build(LineCableSystem,[design,design],[(0.,-1.),(.7,-3.)];
        connections=[Dict("core"=>1),Dict("core"=>2)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,3.),
        frequencies=[50.,1e5],temperature=20.)
    result=compute(problem,Formulation(earth_impedance=:Papadopoulos2010b,
        earth_admittance=:Papadopoulos2010b,
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false));
        options=(trace=true,))
    trace=details(result).trace
    for k in 1:2
        s=2π*problem.frequencies[k]*im
        mu=[μ0,3μ0];epsilon=[ε0,10ε0]
        kx=abs(s)*sqrt(mu[2]*epsilon[2])
        pair=E.EarthPair(1,2,(-1.,-3.),.7,(2,2))
        expected=source_integrals(s,mu,epsilon,[Inf,100.],kx,pair)
        @test trace.Zg[1,2,k] ≈ expected.Z rtol=8e-7
        @test trace.Pg[1,2,k] ≈ expected.P rtol=8e-7
        @test result.Y[:,:,k] ≈ s*inv(trace.P[:,:,k]) rtol=1e-10
        @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k])
    end
end
