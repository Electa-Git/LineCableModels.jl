@testitem "Engine / literature assimilation / Papadopoulos layered F and G witnesses" begin
    using LineCableModels,QuadGK
    E=LineCableModels.Engine;EI=E.EarthImpedance;EA=E.EarthAdmittance
    mu0=4π*1e-7;eps0=8.8541878128e-12
    root(z,s)=iszero(imag(z))&&real(z)<0 ?
        complex(0.,sign(imag(s))*sqrt(-real(z))) : sqrt(z)
    function printed(year,s,rho,epsilon,mu,kx,d,pair)
        g=s.*mu.*(inv.(rho).+s.*epsilon)
        hi,hj=abs.(pair.heights);x=pair.separation
        function kernels(u)
            a=root.(u^2 .+ g .+ kx^2,Ref(s))
            at(i)=a[i+1];mt(i)=mu[i+1];gt(i)=g[i+1]
            S(m,n)=mt(n)*at(m)+mt(m)*at(n)
            D(m,n)=mt(m)*at(n)-mt(n)*at(m)
            A(m,n)=at(n)*gt(m)*mt(n)+at(m)*gt(n)*mt(m)
            Delta(m,n)=at(n)*gt(m)*mt(n)-at(m)*gt(n)*mt(m)
            q=exp(-2at(1)*d)
            if year==2009
                # The 2009 lowercase d has the opposite sign to 2011 D.
                sm=S(0,1)*S(1,2)+D(0,1)*D(1,2)*q
                se=A(0,1)*A(1,2)+Delta(0,1)*Delta(1,2)*q
                F=mt(1)*(S(1,2)-D(1,2)*q)/sm
                numerator=mt(0)*mt(1)*(gt(0)-gt(1))*
                    (S(1,2)-D(1,2)*q)*(A(1,2)-Delta(1,2)*q)-
                    4mt(0)*mt(1)^2*mt(2)*at(1)^2*gt(0)*(gt(2)-gt(1))*q
                G=u*numerator/(sm*se)
                return [F,G]*exp(-u*(hi+hj))*cos(x*u)
            end
            sm=S(1,0)*S(2,1)+D(1,0)*D(2,1)*q
            se=A(1,0)*A(1,2)-Delta(1,0)*Delta(1,2)*q
            e(z)=exp(-at(1)*z)
            F=(S(1,0)*S(2,1)*e(abs(hi-hj))+
                S(1,0)*D(2,1)*e(2d-hi-hj)-
                D(1,0)*S(2,1)*e(hi+hj)-
                D(1,0)*D(2,1)*e(2d-abs(hi-hj)))/(at(1)*sm)
            G1=mt(1)*mt(2)*(gt(1)-gt(2))*
                (S(1,0)*A(1,0)*e(2d-hi-hj)-D(1,0)*A(1,0)*e(2d+hi-hj))
            G2=mt(1)*mt(2)*(gt(1)-gt(2))*
                (S(1,0)*Delta(1,0)*e(2d+hj-hi)-D(1,0)*Delta(1,0)*e(2d+hi+hj))
            G3=mt(1)*mt(0)*(gt(1)-gt(0))*
                (S(2,1)*Delta(1,2)*e(2d+hi-hj)+D(2,1)*Delta(1,2)*e(4d-hi-hj))
            G4=mt(1)*mt(0)*(gt(1)-gt(0))*
                (S(2,1)*A(1,2)*e(hi+hj)+D(2,1)*A(1,2)*e(2d+hj-hi))
            return [F,2at(1)*(G1+G2+G3+G4)/(sm*se)]*cos(x*u)
        end
        cut=g[1]+kx^2
        point=iszero(imag(cut))&&real(cut)<0 ? sqrt(-real(cut)) : 0.
        bounds=point>0 ? (0.,point,Inf) : (0.,Inf)
        integral=quadgk(kernels,bounds...;rtol=2e-9)[1]
        if year==2009
            ideal=log(hypot(x,hi+hj)/hypot(x,hi-hj))
            return (Z=s*mu[1]/(2π)*(ideal+2integral[1]),
                P=(ideal+2sum(integral))/(2π*epsilon[1]))
        end
        return (Z=s*mu[2]/(2π)*integral[1],
            P=s/(2π*(1/rho[2]+s*epsilon[2]))*sum(integral))
    end
    for year in (2009,2011),frequency in (50.,1e5,1e7),permeable in (false,true)
        rho=[Inf,100.,500.];epsilon=eps0*[1.,10.,4.]
        mu=mu0*(permeable ? [1.,3.,.7] : ones(3))
        s=2π*frequency*im
        reference=year==2009 ? 1 : 2
        kx=abs(s)*sqrt(mu[reference]*epsilon[reference])
        id=Symbol("Papadopoulos",year)
        positions=year==2009 ? (10.,15.) : (-.5,-1.5)
        layers=year==2009 ? (1,1) : (2,2)
        for self in (true,false)
            pair=E.EarthPair(1,self ? 1 : 2,
                self ? (positions[1],positions[1]) : positions,
                self ? .02 : .7,layers)
            # Unequal depths ensure convergence of the literal F integral;
            # self gets an independent small-radius direct-Bessel witness
            # through the already checked homogeneous/stratified reduction.
            expected=self ? nothing : printed(year,s,rho,epsilon,mu,kx,3.,pair)
            for (family,field) in ((EI,:Z),(EA,:P))
                leaf=family.Formula(id)(rho,epsilon,mu,s,nothing,nothing,[Inf,3.,Inf])
                kind=self ? Val(:self) : Val(:mutual)
                value=leaf(kind,pair)
                @test isfinite(value)
                self || @test value ≈ getproperty(expected,field) rtol=2e-6
                reversed=E.EarthPair(pair.column,pair.row,reverse(pair.heights),
                    pair.separation,pair.layers)
                @test leaf(kind,reversed) ≈ value rtol=3e-8
                negative=family.Formula(id)(rho,epsilon,mu,-s,nothing,nothing,[Inf,3.,Inf])
                @test negative(kind,pair) ≈ conj(value) rtol=3e-8
                explicit=family.Formula(id)(rho,epsilon,mu,s,complex(kx),nothing,[Inf,3.,Inf])
                @test explicit(kind,pair) ≈ value rtol=3e-8
            end
        end
    end
    for T in (Float32,Float64,BigFloat),family in (EI,EA),year in (2009,2011)
        rho=T[Inf,100,500];epsilon=T(eps0)*T[1,10,4];mu=T(mu0)*T[1,3,1]
        s=complex(zero(T),100T(π))
        heights=year==2009 ? (T(10),T(15)) : (T(-.5),T(-1.5))
        layers=year==2009 ? (1,1) : (2,2)
        pair=E.EarthPair(1,2,heights,T(.7),layers)
        leaf=family.Formula(Symbol("Papadopoulos",year))(
            rho,epsilon,mu,s,nothing,nothing,T[Inf,3,Inf])
        value=leaf(Val(:mutual),pair)
        @test value isa Complex{T}
        @test isfinite(value)
        if year==2011
            @test_throws DomainError leaf(Val(:mutual),
                E.EarthPair(1,2,(T(-.5),T(-4)),T(.7),(2,2)))
        end
    end
end
