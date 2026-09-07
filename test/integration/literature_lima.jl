@testitem "Engine / literature assimilation / Lima closed and asymptotic forms" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0,ε0=4π*1e-7,8.8541878128e-12
    @test formula_id(EI.Formula(:Theodoulidis2015))===:Lima2012
    for frequency in (1.0,50.0,1e4,1e6), resistivity in (10.0,1000.0)
        s=complex(0.0,2π*frequency)
        γ=sqrt(s*μ0*(inv(resistivity)+s*10ε0))
        f=EI.Formula(:Lima2012)([Inf,resistivity],[ε0,10ε0],[μ0,μ0],s,nothing)
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(10.0,10.0),0.01,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(10.0,10.0),3.0,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(10.0,15.0),30.0,(1,1))))
            hi,hj=pair.heights; H=hi+hj
            x=kind===Val(:self) ? 0.0 : pair.separation
            ideal=kind===Val(:self) ? log(H/pair.separation) :
                log(hypot(H,x)/hypot(hi-hj,x))
            integral=quadgk(t->exp(-H*t)*cos(x*t)/(t+sqrt(t^2+γ^2)),
                0.0,Inf;rtol=1e-11)[1]
            expected=s*μ0/(2π)*(ideal+2integral)
            @test f(kind,pair) ≈ expected rtol=1e-8
            reversed=E.EarthPair(pair.column,pair.row,reverse(pair.heights),
                pair.separation,reverse(pair.layers))
            @test f(kind,reversed) ≈ expected rtol=1e-8
        end
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(-1.0,-1.0),0.02,(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-1.0,-2.0),0.5,(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-1.0,-2.0),0.0,(2,2))))
            hi,hj=abs.(pair.heights); H=hi+hj; x=pair.separation
            d=hypot(x,hi-hj); D=hypot(x,H)
            expected=s*μ0/(2π)*(besselk(0,γ*d)+(H^2-x^2)/D^2*besselk(2,γ*D)-
                2*(H^2-x^2)/(γ^2*D^4)*(1+H*γ)*exp(-H*γ))
            @test f(kind,pair) ≈ expected rtol=3e-8
            @test isfinite(f(kind,pair))
        end
        @test_throws DomainError f(Val(:mutual),E.EarthPair(1,2,(-1.0,-1.0),2.0,(2,2)))
        @test_throws ArgumentError f(Val(:mutual),E.EarthPair(1,2,(10.0,-1.0),1.0,(1,2)))
        legacy=EI.Formula(:Theodoulidis2015)([Inf,resistivity],[ε0,10ε0],[μ0,μ0],s,nothing)
        reference=EI.Formula(:Lima2012;displacement_current=false)(
            [Inf,resistivity],[ε0,10ε0],[μ0,μ0],s,nothing)
        pair=E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))
        @test legacy(Val(:mutual),pair) == reference(Val(:mutual),pair)
    end
    # The independent integral is the parent of the exact Struve identity.
    for z in (0.2+0.1im,1.0+1.0im,3.0+4.0im,0.2+8.0im,-1.0+2.0im)
        actual=EI.closed_form_term(Val(:Lima2012),z)
        # Choose gamma and a positive-height complex distance with gamma*distance=z.
        γ=sqrt(im); distance=z/γ
        real(distance)>0 || continue
        source=quadgk(t->exp(-distance*t)/(t+sqrt(t^2+γ^2)),0.0,Inf;rtol=1e-11)[1]
        @test actual ≈ source rtol=1e-8
    end
    for T in (Float32,Float64,BigFloat)
        s=complex(zero(T),T(100)*T(π))
        f=EI.Formula(:Lima2012)(T[Inf,100],T[ε0,10ε0],T[μ0,μ0],s,nothing)
        for (heights,layers) in (((T(10),T(15)),(1,1)),((T(-1),T(-2)),(2,2)))
            pair=E.EarthPair(1,2,heights,T(0.5),layers)
            @test f(Val(:mutual),pair) isa Complex{T}
        end
    end
end
