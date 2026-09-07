@testitem "Engine / literature assimilation / Sunde layered boundaries" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    @test EI.formula_id(EI.Formula(:Sunde1968)) === :Sunde1949
    for frequency in (50.0,10000.0), depth in (2.0,100.0)
        s=complex(0.0,2π*frequency)
        args=([Inf,100.0,500.0],[ε0,10ε0,8ε0],fill(μ0,3),s,nothing,nothing,[Inf,depth,Inf])
        approx=EI.Formula(:Sunde1949;approximation=:large_spacing)(args...)
        exact=EI.Formula(:Sunde1949)(args...)
        g1=sqrt(s*μ0/100); g2=sqrt(s*μ0/500)
        e=exp(-2g1*depth)
        effective=g1*(g1+g2-(g1-g2)*e)/(g1+g2+(g1-g2)*e)
        for (H,x) in ((1e5,0.0),(1e5,2e4),(2e5,1e5))
            pair=E.EarthPair(1,2,(H/3,2H/3),x,(1,1))
            radius=hypot(H,x); theta=atan(x,H)
            correction=s*μ0/π*(cos(theta)/(effective*radius)-
                cos(2theta)/(effective*radius)^2)
            ideal=s*μ0/(2π)*log(radius/hypot(H/3,x))
            @test approx(Val(:mutual),pair) ≈ ideal+correction rtol=1e-12
            @test approx(Val(:mutual),pair) ≈ exact(Val(:mutual),pair) rtol=3e-5
        end
        self=E.EarthPair(1,1,(10.0,10.0),0.02,(1,1))
        @test approx(Val(:self),self) ≈ exact(Val(:self),self) rtol=1e-12
    end
    function model(rho,depths,s; epsilon=zeros(length(rho)))
        return EI.Formula(:Sunde1949; displacement_current=true)(
            [Inf;rho],[ε0;epsilon],fill(μ0,length(rho)+1),s,
            nothing,nothing,[Inf;depths;Inf]
        )
    end
    function explicit_two(λ,g1,g2,d)
        a=sqrt(λ^2+g1); b=sqrt(λ^2+g2); e=exp(-2a*d)
        return (a+b+(a-b)*e)/((a+b)*(λ+a)+(a-b)*(λ-a)*e)
    end
    for frequency in (1.0,50.0,10000.0,1e6),
        (rho1,rho2,depth) in ((100.0,500.0,2.0),(1000.0,10.0,50.0),(100.0,Inf,5.0))
        s=complex(0.0,2π*frequency)
        two=model([rho1,rho2],[depth],s)
        split=model([rho1,rho1,rho2],[depth/3,2depth/3],s)
        for λ in (0.0,1e-5,0.01,1.0,10.0)
            literal=explicit_two(λ,s*μ0/rho1,s*μ0/rho2,depth)
            @test EI.spectral_kernel(Val(:Sunde1949),λ,two.state) ≈ literal rtol=1e-11
            @test EI.spectral_kernel(Val(:Sunde1949),λ,split.state) ≈ literal rtol=1e-11
        end
        for pair in (E.EarthPair(1,1,(10.0,10.0),0.02,(1,1)),
                E.EarthPair(1,2,(10.0,15.0),3.0,(1,1)),
                E.EarthPair(1,2,(10.0,15.0),0.0,(1,1)))
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            @test split(kind,pair) ≈ two(kind,pair) rtol=3e-8
        end
        # Iwamoto August 1958 appendix (付1), including its printed MKS
        # prefactor without an extra j. The integral is the ground correction.
        scale=sqrt(imag(s)*μ0/rho1); h=10.0
        normalized=quadgk(0.0,Inf;rtol=1e-10) do t
            a=sqrt(t^2+im); b=sqrt(t^2+im*rho1/rho2)
            r12=(a-b)/(a+b); r10=(a-t)/(a+t)
            e=exp(-2scale*depth*a)
            im/(a+t)*(1+r12*e)/(1-r12*r10*e)*exp(-2scale*h*t)
        end[1]
        expected=imag(s)*μ0/π*normalized
        pair=E.EarthPair(1,1,(h,h),0.02,(1,1))
        correction=two(Val(:self),pair)-s*μ0/(2π)*log(2h/pair.separation)
        @test correction ≈ expected rtol=3e-8
    end
    for count in 2:5, frequency in (1.0,50.0,10000.0,1e6)
        s=complex(0.0,2π*frequency)
        homogeneous=model([100.0],Float64[],s;epsilon=[10ε0])
        layered=model(fill(100.0,count),fill(3.0,count-1),s;epsilon=fill(10ε0,count))
        for pair in (E.EarthPair(1,1,(10.0,10.0),0.02,(1,1)),
                E.EarthPair(1,2,(10.0,15.0),3.0,(1,1)))
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            @test layered(kind,pair) ≈ homogeneous(kind,pair) rtol=3e-8
        end
    end
end
