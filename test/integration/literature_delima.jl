@testitem "Engine / literature assimilation / De Lima complex soil" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance; FD=LineCableModels.EarthProps.FD
    μ0,ε0=4π*1e-7,8.8541878128e-12
    Δ,α=8.92028e-3,0.71603
    β=Δ*cot(π*α/2)/(1e-6*(2π*1e6)^α)
    law=FD.Formula(:Portela1999;beta=β,exponent=α)
    initial=LineCableModels.EarthProps.EarthMaterial(inv(84.16e-6),1.0,1.0)
    for frequency in (1.0,50.0,1e4,1e6)
        s=complex(0.0,2π*frequency)
        material=constitutive(law,initial,frequency)
        κ=84.16e-6+Δ*(frequency/1e6)^α*(cot(α*π/2)+im)
        @test inv(material.rho)+s*ε0*material.eps_r ≈ κ rtol=1e-12
        rho=[Inf,material.rho]; ε=[ε0,ε0*material.eps_r]; μ=[μ0,μ0]
        γ=sqrt(s*μ0*κ)
        functor=EI.Formula(:DeLima2007)(rho,ε,μ,s,nothing)
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(10.0,10.0),0.01,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))),
            (Val(:self),E.EarthPair(1,1,(-1.0,-1.0),0.02,(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-1.0,-2.0),0.5,(2,2))))
            hi,hj=abs.(pair.heights); x=pair.separation; H=hi+hj
            if pair.layers[1]==1
                lateral=kind===Val(:self) ? 0.0 : x
                correction=quadgk(t->exp(-H*t)*cos(lateral*t)/
                    (t+sqrt(t^2+γ^2)),0.0,Inf;rtol=1e-10)[1]
                base=kind===Val(:self) ? log(H/x) :
                    log(hypot(x,H)/hypot(x,hi-hj))
                expected=s*μ0/(2π)*(base+2correction)
            else
                correction=quadgk(t->exp(-H*sqrt(t^2+γ^2))*cos(x*t)/
                    (t+sqrt(t^2+γ^2)),0.0,Inf;rtol=1e-10)[1]
                # Appendix (30), not the inconsistent main-text K1 image.
                expected=s*μ0/(2π)*(besselk(0,γ*hypot(x,hi-hj))-
                    besselk(0,γ*hypot(x,H))+2correction)
            end
            @test functor(kind,pair) ≈ expected rtol=3e-8
            @test isfinite(functor(kind,pair))
        end
        @test_throws ArgumentError functor(Val(:mutual),
            E.EarthPair(1,2,(10.0,-1.0),1.0,(1,2)))
    end
    for T in (Float32,Float64,BigFloat)
        s=complex(zero(T),T(100)*T(π))
        rho=T[Inf,100]; ε=T[ε0,0]; μ=T[μ0,μ0]
        f=EI.Formula(:DeLima2007)(rho,ε,μ,s,nothing)
        pair=E.EarthPair(1,2,(T(-1),T(-2)),T(0.5),(2,2))
        parent=EI.Formula(:Pollaczek1926)(rho,ε,μ,s,nothing)
        @test f(Val(:mutual),pair) ≈ parent(Val(:mutual),pair)
        @test f(Val(:mutual),pair) isa Complex{T}
    end
end
