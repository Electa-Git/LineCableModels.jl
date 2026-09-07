@testitem "Engine / literature assimilation / Xue exact four-layer coefficients" begin
    using LineCableModels, QuadGK, LinearAlgebra
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance; EP=LineCableModels.EarthProps
    mu0=4π*1e-7; ep0=8.8541878128e-12
    # Direct transcription of final publication equations (10)–(44).
    # Indices in this witness start with air=1; paper interface depths
    # are positive cumulative distances below the interface.
    function xue_appendix(lambda,g,mu,thickness)
        a=sqrt.(lambda^2 .+ g .- g[1])
        d=cumsum(thickness)
        function factors(n,m)
            ep=exp((-a[n]+a[m])*d[n-1]); eq=exp((a[n]+a[m])*d[n-1])
            ew=inv(eq); ev=inv(ep)
            ratio=g[n]/g[m]; magnetic=mu[m]*a[n]/(mu[n]*a[m])
            P=ratio/2*(1+magnetic)*ep; Q=ratio/2*(1-magnetic)*eq
            W=ratio/2*(1-magnetic)*ew; V=ratio/2*(1+magnetic)*ev
            PG=(a[n]/a[m]+mu[m]/mu[n]*ratio)/2*ep
            QG=(a[n]/a[m]-mu[m]/mu[n]*ratio)/2*eq
            WG=(a[n]/a[m]-mu[m]/mu[n]*ratio)/2*ew
            VG=(a[n]/a[m]+mu[m]/mu[n]*ratio)/2*ev
            return (;P,Q,W,V,PG,QG,WG,VG,ep,eq,ew,ev)
        end
        f=factors(2,3); t=factors(3,4)
        S=g[1]/g[2]*(1+mu[2]*lambda/(mu[1]*a[2]))
        Sp=g[1]/g[2]*(1-mu[2]*lambda/(mu[1]*a[2]))
        D1=(f.P*S+f.Q*Sp)/2; D2=lambda*(f.P*Sp+f.Q*S)/2
        D3=(f.W*S+f.V*Sp)/2; D4=lambda*(f.W*Sp+f.V*S)/2
        T1=t.P*D1+t.Q*D3; T2=t.P*D2+t.Q*D4
        T3=t.W*D1+t.V*D3; T4=t.W*D2+t.V*D4
        T5=T2/lambda; T6=T4/lambda
        A1=(1-mu[5]*a[4]/(mu[4]*a[5]))*exp(-a[4]*d[3])
        A2=(1+mu[5]*a[4]/(mu[4]*a[5]))*exp(a[4]*d[3])
        F=(A1*(T5-T1)+A2*(T6-T3))/(A1*T2+A2*T4)
        A3=(g[1]/g[2]-1)/a[2]
        B=g[1]*mu[2]/(g[2]*mu[1])+lambda/a[2]
        Bp=g[1]*mu[2]/(g[2]*mu[1])-lambda/a[2]
        exp2=exp(2a[3]*d[1])
        D5=A3/2*(f.PG+f.QG)+S/(4a[3])*(f.P+f.W*exp2-f.ep)+
            Sp/(4a[3])*(f.Q+f.V*exp2-f.eq)
        D6=A3/2*(f.PG+f.QG)+S/(4a[3])*(f.Q+f.V*exp2-f.eq)+
            Sp/(4a[3])*(f.P+f.W*exp2-f.ep)
        D7=(f.PG*Bp-f.QG*B)/2
        D8=-A3/2*(f.WG+f.VG)-S/(4a[3])*(f.W+f.P/exp2-f.ew)-
            Sp/(4a[3])*(f.V+f.Q/exp2-f.ev)
        D9=-A3/2*(f.VG+f.WG)-S/(4a[3])*(f.V+f.Q/exp2-f.ev)-
            Sp/(4a[3])*(f.W+f.P/exp2-f.ew)
        D10=(f.VG*B-f.WG*Bp)/2
        exp3=exp(2a[4]*d[2])
        D11=t.VG*D8-t.WG*D5+t.ew/(2a[4])*D1+t.ev/(2a[4])*D3-
            T1/(2a[4])/exp3-T3/(2a[4])
        D12=t.VG*D9-t.WG*D6+t.ew/(2a[4]*lambda)*D2+t.ev/(2a[4]*lambda)*D4-
            T5/(2a[4])/exp3-T6/(2a[4])
        D13=t.VG*D10-t.WG*D7
        D14=t.PG*D5-t.QG*D8-t.ep/(2a[4])*D1-t.eq/(2a[4])*D3+
            T3/(2a[4])*exp3+T1/(2a[4])
        D15=t.PG*D6-t.QG*D9-t.ep/(2a[4]*lambda)*D2-t.eq/(2a[4]*lambda)*D4+
            T6/(2a[4])*exp3+T5/(2a[4])
        D16=t.PG*D7-t.QG*D10
        ratio=g[4]/g[5]; magnetic=ratio*mu[5]/mu[4]*a[5]
        low=exp(-a[4]*d[3]); high=inv(low)
        Tp1=(T1-a[4]*D14-ratio*T1+magnetic*D14)*low+
            (T3+a[4]*D11-ratio*T3+magnetic*D11)*high
        Tp2=(T5-a[4]*D15-ratio*T5+magnetic*D15)*low+
            (T6+a[4]*D12-ratio*T6+magnetic*D12)*high
        Tp3=(a[4]*D16-magnetic*D16)*low-(a[4]*D13+magnetic*D13)*high
        Tp4=-(A1*T1+A2*T3)/(A1*T5+A2*T6)
        G=F+(Tp1+Tp2*Tp4)/Tp3
        return (;F,G)
    end

    function original_integrals(s,rho,epsilon,mu,thickness,H,x,d)
        g=Complex{BigFloat}.(s.*mu.*(inv.(rho).+s.*epsilon))
        mb=BigFloat.(mu); hb=BigFloat.(thickness[2:4])
        # Source appendix evaluated at high precision to retain F+electric
        # cancellation. The finite upper bound gives exp(-80) attenuation.
        integral=quadgk(0.,80/H;rtol=1e-10) do lambda
            source=xue_appendix(BigFloat(lambda),g,mb,hb)
            ComplexF64[source.F,source.G]*exp(-lambda*H)*cos(lambda*x)
        end[1]
        ideal=log(hypot(H,x)/d)
        return (Z=s*mu[1]/(2π)*(ideal+integral[1]),
            P=(ideal+integral[2])/(2π*epsilon[1]))
    end
    for f in (1.,50.,1e5,1e7),lambda in (1e-5,.001,.1,2.),magnetic in (false,true)
        s=2π*f*im; rho=[Inf,100.,500.,50.,200.]; epsilon=ep0.*[1.,10.,6.,20.,8.]
        mu=mu0.*(magnetic ? [1.,2.,3.,4.,5.] : ones(5)); thickness=[Inf,2.,3.,5.,Inf]
        functor=EA.Formula(:Xue2021;evaluation=:exact)(rho,epsilon,mu,s,nothing,nothing,thickness)
        g=Complex{BigFloat}.(s.*mu.*(inv.(rho).+s.*epsilon))
        reference=xue_appendix(BigFloat(lambda),g,BigFloat.(mu),BigFloat.(thickness[2:4]))
        actual=EA.earth_potential_coefficient(Val(:Xue2021),Val(:kernel),lambda,functor.state)
        @test actual.F ≈ reference.F rtol=1e-10
        @test actual.G ≈ reference.G rtol=2e-9
        magnetic_functor=EI.Formula(:Xue2021)(rho,epsilon,mu,s,nothing,nothing,thickness)
        @test 2EI._layered_overhead_kernel(lambda,magnetic_functor.state) ≈ reference.F rtol=1e-10
    end
    for f in (50.,1e5),magnetic in (false,true)
        s=2π*f*im; rho=[Inf,100.,500.,50.,200.]; epsilon=ep0.*[1.,10.,6.,20.,8.]
        mu=mu0.*(magnetic ? [1.,2.,3.,4.,5.] : ones(5)); thickness=[Inf,2.,3.,5.,Inf]
        Z=EI.Formula(:Xue2021)(rho,epsilon,mu,s,nothing,nothing,thickness)
        P=EA.Formula(:Xue2021;evaluation=:exact)(rho,epsilon,mu,s,nothing,nothing,thickness)
        for pair in (E.EarthPair(1,1,(10.,10.),.01,(1,1)),E.EarthPair(1,2,(10.,15.),3.,(1,1)))
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            H=sum(pair.heights); x=kind===Val(:self) ? 0. : pair.separation
            d=kind===Val(:self) ? pair.separation : hypot(x,pair.heights[1]-pair.heights[2])
            reference=original_integrals(s,rho,epsilon,mu,thickness,H,x,d)
            @test Z(kind,pair) ≈ reference.Z rtol=3e-8
            @test P(kind,pair) ≈ reference.P rtol=3e-8
            reversed=E.EarthPair(pair.column,pair.row,reverse(pair.heights),pair.separation,pair.layers)
            @test Z(kind,reversed) ≈ Z(kind,pair) rtol=1e-12
            @test P(kind,reversed) ≈ P(kind,pair) rtol=1e-12
            @test EI.Formula(:Xue2021)(rho,epsilon,mu,-s,nothing,nothing,thickness)(kind,pair) ≈ conj(Z(kind,pair))
            @test EA.Formula(:Xue2021;evaluation=:exact)(rho,epsilon,mu,-s,nothing,nothing,thickness)(kind,pair) ≈ conj(P(kind,pair))
        end
    end
    for T in (Float32,Float64,BigFloat)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),T(100)*T(π)); rho=T[Inf,100,500,50,200]
        epsvec=epsilon.*T[1,10,6,20,8]; muvec=mu.*T[1,2,3,4,5]; h=T[Inf,2,3,5,Inf]
        pair=E.EarthPair(1,2,(T(10),T(15)),T(3),(1,1))
        for recipe in (EI.Formula(:Xue2021),EA.Formula(:Xue2021;evaluation=:exact))
            value=recipe(rho,epsvec,muvec,s,nothing,nothing,h)(Val(:mutual),pair)
            @test value isa Complex{T}
            @test isfinite(value)
        end
    end
    s=100π*im; pair=E.EarthPair(1,2,(10.,15.),3.,(1,1))
    rho=[Inf,100.,100.,100.,100.]; epsilon=ep0.*[1.,10.,10.,10.,10.]; mu=fill(mu0,5)
    for h in ([Inf,0.,0.,0.,Inf],[Inf,2.,3.,5.,Inf],[Inf,1e6,1e6,1e6,Inf])
        z=EI.Formula(:Xue2021)(rho,epsilon,mu,s,nothing,nothing,h)(Val(:mutual),pair)
        p=EA.Formula(:Xue2021;evaluation=:exact)(rho,epsilon,mu,s,nothing,nothing,h)(Val(:mutual),pair)
        @test z ≈ EI.Formula(:Wise1934)(rho[1:2],epsilon[1:2],mu[1:2],s,nothing)(Val(:mutual),pair) rtol=3e-8
        @test p ≈ EA.Formula(:Wise1948)(rho[1:2],epsilon[1:2],mu[1:2],s,nothing)(Val(:mutual),pair) rtol=3e-8
    end
    # Two-layer limit: split the bottom medium into identical layers.
    rho=[Inf,100.,500.,500.,500.]; epsilon=ep0.*[1.,10.,6.,6.,6.]; mu=mu0.*[1.,2.,3.,3.,3.]
    h=[Inf,2.,3.,5.,Inf]
    p=EA.Formula(:Xue2021;evaluation=:exact)(rho,epsilon,mu,s,nothing,nothing,h)(Val(:mutual),pair)
    parent=EA.Formula(:Papadopoulos2009)(rho[1:3],epsilon[1:3],mu[1:3],s,nothing,nothing,[Inf,2.,Inf])
    @test p ≈ parent(Val(:mutual),pair) rtol=3e-8

    for recipe in (EI.Formula(:Xue2021),EA.Formula(:Xue2021;evaluation=:exact))
        @test_throws DimensionMismatch recipe(rho[1:2],epsilon[1:2],mu[1:2],s,nothing)
        @test_throws ArgumentError recipe(rho,epsilon,mu,s,zero(s),nothing,h)
        @test_throws DomainError recipe(rho,epsilon,mu,s,nothing,nothing,[Inf,-1.,2.,3.,Inf])
        @test_throws DomainError recipe(rho,epsilon,mu,zero(s),nothing,nothing,h)
        functor=recipe(rho,epsilon,mu,s,nothing,nothing,h)
        @test_throws ArgumentError functor(Val(:mutual),E.EarthPair(1,2,(10.,-1.),3.,(1,2)))
    end
    @test_throws ArgumentError EA.Formula(:Xue2021;evaluation=:invalid)
    @test E.media(EA.Formula(:Xue2021))===Val(:homogeneous)
    @test E.media(EA.Formula(:Xue2021;evaluation=:exact))===Val(:stratified)
    @test EA.propagation(EA.Formula(:Xue2021))===Val(:zero)
    copper=Material(:conductor,1.7241e-8,1.,1.,20.,0.)
    wire=build(CableDesign,"Xue exact wire",Group(:core,Region(:core,Disk(.01),copper)))
    system=build(LineCableSystem,[wire,wire],[(0.,10.),(3.,15.)];connections=[Dict("core"=>1),Dict("core"=>2)])
    earth=build(EarthModel,(EP.EarthLayer(100.,10.,2.,2.),EP.EarthLayer(500.,6.,3.,3.),
        EP.EarthLayer(50.,20.,4.,5.),EP.EarthLayer(200.,8.,5.)))
    problem=LineParametersProblem(system;earth_props=earth,frequencies=[50.,1e5])
    formulation=Formulation(earth_impedance=:Xue2021,earth_admittance=formula(:Xue2021;evaluation=:exact),
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    result=compute(problem,formulation;options=(trace=true,))
    trace=details(result).trace
    for (k,f) in enumerate(problem.frequencies)
        freq_s=2π*f*im; expectedZ=zeros(ComplexF64,2,2); expectedP=similar(expectedZ)
        for i in 1:2,j in i:2
            H=[10.,15.][i]+[10.,15.][j]; x=i==j ? 0. : 3.; d=i==j ? .01 : hypot(3.,5.)
            ref=original_integrals(freq_s,[Inf,100.,500.,50.,200.],ep0.*[1.,10.,6.,20.,8.],
                mu0.*[1.,2.,3.,4.,5.],[Inf,2.,3.,5.,Inf],H,x,d)
            expectedZ[i,j]=expectedZ[j,i]=ref.Z; expectedP[i,j]=expectedP[j,i]=ref.P
        end
        @test trace.Zg[:,:,k] ≈ expectedZ rtol=3e-8
        @test trace.Pg[:,:,k] ≈ expectedP rtol=3e-8
        @test result.Y[:,:,k] ≈ freq_s*inv(expectedP) rtol=3e-8
        @test minimum(eigvals(Symmetric(real.(result.Y[:,:,k]))))>=-1e-10
    end
end
