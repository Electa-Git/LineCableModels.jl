@testitem "Engine / literature assimilation / Pettersson image prescriptions" begin
    using LineCableModels
    E=LineCableModels.Engine;EI=E.EarthImpedance;EA=E.EarthAdmittance
    mu0=4π*1e-7;eps0=8.8541878128e-12
    # Independent source (10), (11), (14), (15), with medium 1 local
    # to the wire; the engine helper is not used by this witness.
    function printed(s,rho,epsilon,pair)
        self=pair.row==pair.column
        h,y=abs.(pair.heights)
        interface=h==y==0
        local_medium=interface ? 1 : first(pair.layers)
        other=3-local_medium
        k1=1/rho[local_medium]+s*epsilon[local_medium]
        k2=1/rho[other]+s*epsilon[other]
        n2=k2/k1
        beta=sqrt(s*mu0*(k2-k1))
        if interface
            if self
                P=log(2sqrt(im)/(beta*pair.separation))
                dQ=im*sqrt(2)*(n2+1)/(n2+im)/(beta*pair.separation)
                imag(dQ)>0 && (dQ=-dQ)
                Q=2/(n2+1)*log(dQ)
                return (Z=s*mu0/(2π)*P,potential=s/(2π*k1)*Q)
            end
            d=pair.separation;mirror=d;L=0.
            dp=sqrt((sqrt(2)*(1+im)/beta)^2+d^2)
            dq=sqrt((im*sqrt(2)*(n2+1)/(n2+im)/beta)^2+d^2)
        elseif self
            d=pair.separation;mirror=2h;L=log(mirror/d)
            dp=2h+2/beta
            dq=2h+(n2+1)/beta
        else
            d=hypot(pair.separation,y-h)
            mirror=hypot(pair.separation,y+h);L=log(mirror/d)
            dp=sqrt((y+h+2/beta)^2+pair.separation^2)
            dq=sqrt((y+h+(n2+1)/beta)^2+pair.separation^2)
        end
        wanted=local_medium==1 ? -1 : 1
        imag(dq)*wanted<0 && (dq=-dq)
        return (Z=s*mu0/(2π)*(L+log(dp/mirror)),
            potential=s/(2π*k1)*(L+2/(n2+1)*log(dq/mirror)))
    end
    pairs=(E.EarthPair(1,1,(10.,10.),.01,(1,1)),
        E.EarthPair(1,2,(10.,12.),5.,(1,1)),
        E.EarthPair(1,1,(-1.,-1.),.01,(2,2)),
        E.EarthPair(1,2,(-1.,-2.),.5,(2,2)),
        E.EarthPair(1,1,(0.,0.),.01,(1,1)),
        E.EarthPair(1,2,(0.,0.),.5,(1,1)))
    for frequency in (50.,5e3,1e6,5e6),airrho in (Inf,1e8)
        s=2π*frequency*im;rho=[airrho,1000.]
        epsilon=[eps0,10eps0];mu=[mu0,mu0]
        for pair in pairs
            expected=printed(s,rho,epsilon,pair)
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            for (family,field) in ((EI,:Z),(EA,:potential))
                functor=family.Formula(:Pettersson1994)(rho,epsilon,mu,s,nothing)
                value=functor(kind,pair)
                @test value ≈ getproperty(expected,field) rtol=2e-12
                negative=family.Formula(:Pettersson1994)(rho,epsilon,mu,-s,nothing)
                @test negative(kind,pair) ≈ conj(value) rtol=2e-12
                reversed=E.EarthPair(pair.column,pair.row,reverse(pair.heights),
                    pair.separation,reverse(pair.layers))
                @test functor(kind,reversed) ≈ value rtol=2e-12
            end
        end
    end
    for T in (Float32,Float64,BigFloat),family in (EI,EA),pair in pairs
        rho=T[Inf,1000];epsilon=T[eps0,10eps0];mu=T[mu0,mu0]
        s=complex(zero(T),2T(π)*T(5000))
        leaf=family.Formula(:Pettersson1994)(rho,epsilon,mu,s,nothing)
        typed=E.EarthPair(pair.row,pair.column,T.(pair.heights),T(pair.separation),pair.layers)
        kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
        value=leaf(kind,typed)
        @test value isa Complex{T}
        @test isfinite(value)
    end
    for family in (EI,EA)
        leaf=family.Formula(:Pettersson1994)([Inf,1000.],[eps0,10eps0],
            [mu0,mu0],2π*5e3*im,nothing)
        @test_throws ArgumentError leaf(Val(:mutual),E.EarthPair(1,2,(1.,-1.),1.,(1,2)))
        @test_throws ArgumentError leaf(Val(:mutual),E.EarthPair(1,2,(0.,1.),1.,(1,1)))
        @test_throws DomainError leaf(Val(:self),E.EarthPair(1,1,(0.,0.),1e6,(1,1)))
        # Check continuity through the principal-square-root sign change.
        for frequency in (1e3,1e5,1e7)
            values=map((1-1e-7,1+1e-7)) do ratio
                f=family.Formula(:Pettersson1994)([Inf,1000.],[eps0,10eps0],
                    [mu0,mu0],2π*frequency*ratio*im,nothing)
                f(Val(:mutual),pairs[2])
            end
            @test abs(values[2]/values[1]-1)<1e-5
        end
    end
end
