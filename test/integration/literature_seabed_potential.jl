@testitem "Engine / literature assimilation / three-medium seabed potential" begin
    using LineCableModels, LinearAlgebra, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EA=E.EarthAdmittance; EP=LineCableModels.EarthProps
    mu0=4π*1e-7; ep0=8.8541878128e-12
    # Equations (25)–(27), with weighted potentials so continuity rows
    # remain well scaled. The imposed incident x-potential is unity.
    function boundary(a,g,mu,h)
        a0,a1,a2=a; g0,g1,g2=g; m0,m1,m2=mu; e=exp(-a1*h)
        X=ComplexF64[
            1 -1 -e 0
            -a0/m0 -a1/m1 a1*e/m1 0
            0 e 1 -1
            0 a1*e/m1 -a1/m1 -a2/m2
        ]
        x=X\ComplexF64[0,0,1,-a2/m2]
        V=ComplexF64[
            1 -1 -e 0
            -m0*g1*a0 -m1*g0*a1 m1*g0*a1*e 0
            0 e 1 -1
            0 m1*g2*a1*e -m1*g2*a1 -m2*g1*a2
        ]
        rhs=ComplexF64[0,m2*(g1*x[1]-g0*(x[2]+e*x[3])),0,
            m2*(g2*(e*x[2]+x[3])-g1*(1+x[4]))]
        for i in 1:4
            scale=maximum(abs,V[i,:]); V[i,:]/=scale; rhs[i]/=scale
        end
        v=V\rhs
        return x[4],-v[4]
    end
    function source_reference(s,mu,epsilon,rho,hs,di,dj,x)
        g=s.*mu.*(inv.(rho).+s.*epsilon)
        d=hypot(x,di-dj); H=di+dj
        correction=quadgk(0.,Inf;rtol=1e-10) do t
            lambda=t/H
            a=sqrt.(lambda^2 .+g)
            R,G=boundary(a,g,mu,hs)
            (R/a[3]+G)*exp(-a[3]*H)*cos(lambda*x)/H
        end[1]
        return s/(2π*(inv(rho[3])+s*epsilon[3]))*(besselk(0,sqrt(g[3])*d)+correction)
    end
    for f in (1.,50.,1e5),hs in (0.,.1,1.,10.),mur in ((1.,1.,1.),(1.,2.,3.))
        mu=mu0.*collect(mur); epsilon=ep0.*[1.,80.,10.]; rho=[Inf,.2,100.]
        s=2π*f*im
        functor=EA.Formula(:DiLorenzo2023)(rho,epsilon,mu,s,nothing,nothing,[Inf,hs,Inf])
        g=s.*mu.*(inv.(rho).+s.*epsilon)
        for lambda in (.001,.1,10.)
            a=sqrt.(lambda^2 .+g)
            R,G=boundary(a,g,mu,hs)
            kernel=EA.earth_potential_coefficient(Val(:DiLorenzo2023),Val(:kernel),lambda,functor.state)
            @test kernel.R ≈ R rtol=1e-9 atol=1e-12
            @test kernel.G ≈ G rtol=1e-8 atol=1e-12
        end
        for (di,dj,x,kind) in ((.5,1.,.3,Val(:mutual)),(.5,.5,.02,Val(:self)))
            pair=E.EarthPair(1,kind===Val(:self) ? 1 : 2,(-hs-di,-hs-dj),x,(3,3))
            value=functor(kind,pair)
            expected=source_reference(s,mu,epsilon,rho,hs,di,dj,x)
            @test value ≈ expected rtol=3e-8
            negative=EA.Formula(:DiLorenzo2023)(rho,epsilon,mu,-s,nothing,nothing,[Inf,hs,Inf])
            @test negative(kind,pair) ≈ conj(value) rtol=3e-8
            if kind===Val(:mutual)
                reverse=E.EarthPair(2,1,(-hs-dj,-hs-di),x,(3,3))
                @test functor(kind,reverse) ≈ value rtol=1e-12
            end
            if iszero(hs)
                reference=EA.Formula(:MartinsBritto2024)(rho[[1,3]],epsilon[[1,3]],
                    mu[[1,3]],s,0.0im)
                two=E.EarthPair(pair.row,pair.column,(-di,-dj),x,(2,2))
                @test value ≈ reference(kind,two) rtol=3e-8
                if mur==(1.,1.,1.)
                    nonmagnetic=EA.Formula(:Papadopoulos2010b)(rho[[1,3]],epsilon[[1,3]],
                        mu[[1,3]],s,0.0im)
                    @test value ≈ nonmagnetic(kind,two) rtol=3e-8
                end
            end
        end
    end
    for T in (Float32,Float64,BigFloat)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        rho=T[Inf,.2,100]; epsvec=T[epsilon,80epsilon,10epsilon]; muvec=T[mu,mu,mu]
        s=complex(zero(T),T(100)*T(π)); hs=one(T)
        f=EA.Formula(:DiLorenzo2023)(rho,epsvec,muvec,s,nothing,nothing,T[Inf,hs,Inf])
        pair=E.EarthPair(1,2,(T(-1.5),T(-2)),T(.3),(3,3))
        value=f(Val(:mutual),pair)
        reference=source_reference(ComplexF64(s),Float64.(muvec),Float64.(epsvec),
            Float64.(rho),1.,.5,1.,.3)
        @test value isa Complex{T}
        @test value ≈ reference rtol=(T===Float32 ? 1e-4 : 3e-8)
    end
    # Removing either interface must give the same two-medium parent.
    s=100π*im; rho=[Inf,100.,100.]; epsvec=ep0.*[1.,10.,10.]; muvec=fill(mu0,3)
    homogeneous=EA.Formula(:Papadopoulos2010b)(rho[1:2],epsvec[1:2],muvec[1:2],s,0.0im)
    for hs in (.1,1.,10.)
        f=EA.Formula(:DiLorenzo2023)(rho,epsvec,muvec,s,nothing,nothing,[Inf,hs,Inf])
        pair=E.EarthPair(1,2,(-hs-.5,-hs-1.),.3,(3,3))
        two=E.EarthPair(1,2,pair.heights,pair.separation,(2,2))
        @test f(Val(:mutual),pair) ≈ homogeneous(Val(:mutual),two) rtol=3e-8
    end
    f=EA.Formula(:DiLorenzo2023)(rho,epsvec,muvec,s,nothing,nothing,[Inf,1.,Inf])
    @test_throws ArgumentError f(Val(:mutual),E.EarthPair(1,2,(-.5,-2.),.3,(2,3)))
    @test_throws DomainError f(Val(:mutual),E.EarthPair(1,2,(-.5,-2.),.3,(3,3)))
    @test_throws ArgumentError EA.Formula(:DiLorenzo2023)(rho,epsvec,muvec,s,1.0+0im,
        nothing,[Inf,1.,Inf])
    @test_throws DomainError EA.Formula(:DiLorenzo2023)(rho,epsvec,muvec,s,nothing,
        nothing,[Inf,-1.,Inf])

    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    insulation=Material(:insulator,Inf,2.3,1.0,20.0,0.0)
    design=build(CableDesign,"seabed-coax",Stack(
        Group(:core,Region(:core,Disk(.01),copper)),Region(:inner,Shell(.004),insulation),
        Group(:sheath,Region(:sheath,Shell(.001),copper)),Region(:outer,Shell(.002),insulation)))
    system=build(LineCableSystem,[design,design],[(0.,-1.5),(1.,-2.)];
        connections=[Dict("core"=>1,"sheath"=>2),Dict("core"=>3,"sheath"=>4)])
    earth=build(EarthModel,(EP.EarthLayer(.2,80.,1.,1.),EP.EarthLayer(100.,10.,1.)))
    problem=LineParametersProblem(system;earth_props=earth,frequencies=[50.,1e4],temperature=20.)
    formulation=Formulation(earth_impedance=:Tsiamitros2008,earth_admittance=:DiLorenzo2023,
        insulation_admittance=:Ametani1980,
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    result=compute(problem,formulation;options=(trace=true,))
    trace=details(result).trace
    incidence=[1. 0 0 0;-1 1 0 0;0 0 1 0;0 0 -1 1]
    pcs=log(.014/.01)/(2π*2.3ep0); pse=log(.017/.015)/(2π*2.3ep0)
    @test size(result.Y)==(4,4,2)
    @test all(isfinite,result.Z)
    for k in 1:2
        frequency_s=2π*problem.frequencies[k]*im
        external=[source_reference(frequency_s,fill(mu0,3),ep0.*[1.,80.,10.],[Inf,.2,100.],
            1.,i==1 ? .5 : 1.,j==1 ? .5 : 1.,i==j ? .017 : 1.) for i in 1:2,j in 1:2]
        loops=zeros(ComplexF64,4,4)
        loops[1,1]=pcs;loops[3,3]=pcs
        loops[2,2]=pse+external[1,1];loops[4,4]=pse+external[2,2]
        loops[2,4]=external[1,2];loops[4,2]=external[2,1]
        phase=inv(transpose(incidence))*loops*inv(incidence)
        @test trace.Pg[:,:,k] ≈ external rtol=3e-8
        @test trace.P[:,:,k] ≈ phase rtol=3e-8
        @test result.Y[:,:,k] ≈ frequency_s*inv(phase) rtol=3e-8
        @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k]) rtol=1e-12
        @test minimum(eigvals(Symmetric(real.(result.Y[:,:,k]))))>=-1e-10
    end
end
