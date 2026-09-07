@testitem "Engine / literature assimilation / Carson and Iwamoto normalization" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    for frequency in (1.0,50.0,10000.0,1e6), resistivity in (10.0,1000.0)
        ω=2π*frequency; s=complex(0.0,ω); scale=sqrt(ω*μ0/resistivity)
        f=EI.Formula(:Carson1926)([Inf,resistivity],[ε0,10ε0],[μ0,μ0],s,nothing)
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(1.0,1.0),0.005,(1,1))),
            (Val(:self),E.EarthPair(1,1,(1.0,1.0),0.05,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(1.0,2.0),0.0,(1,1))))
            hi,hj=pair.heights; H=hi+hj
            x=kind===Val(:self) ? 0.0 : pair.separation
            # Iwamoto's printed SI/MKS parent (2),(3); rationalize only
            # sqrt(t²+j)-t to avoid subtracting nearly equal real parts.
            J=quadgk(t->im/(sqrt(t^2+im)+t)*exp(-scale*H*t)*cos(scale*x*t),
                0.0,Inf;rtol=1e-11)[1]
            ideal=kind===Val(:self) ? log(H/pair.separation) :
                log(hypot(H,x)/hypot(hi-hj,x))
            expected=ω*1e-7*(2im*ideal+4J)
            @test f(kind,pair) ≈ expected rtol=3e-8
            @test real(f(kind,pair))>=0
        end
        first=E.EarthPair(1,1,(1.0,1.0),0.005,(1,1))
        second=E.EarthPair(1,1,(1.0,1.0),0.05,(1,1))
        correction1=f(Val(:self),first)-s*μ0/(2π)*log(2/first.separation)
        correction2=f(Val(:self),second)-s*μ0/(2π)*log(2/second.separation)
        @test correction1 ≈ correction2 rtol=1e-11
    end
end
