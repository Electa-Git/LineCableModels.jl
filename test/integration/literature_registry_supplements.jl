@testitem "Engine / literature assimilation / registry supplements and legacy constructors" begin
    using LineCableModels, LinearAlgebra
    E=LineCableModels.Engine; EI=E.EarthImpedance
    EA=E.EarthAdmittance; II=E.InternalImpedance
    for T in (Float32,Float64,BigFloat)
        tol=T===Float32 ? T(4e-5) : T(3e-12)
        μ=T(4)*T(π)*T(10)^(-7); ε=T(88541878128)*T(10)^(-22)
        for frequency in T.((0,50,1e5,1e10)), mr in T.((1,10)),
            (a,b) in ((T(.01),T(.012)),(T(.1),T(.101)))
            ρ=T(1.724e-8); s=complex(zero(T),2T(π)*frequency)
            m=sqrt(s*μ*mr/ρ); t=b-a
            # Independent evaluation of the printed hyperbolic expressions.
            wide=Complex{BigFloat}(m); width=BigFloat(t)
            mc=iszero(m) ? inv(width) : wide*coth(wide*width)
            ms=iszero(m) ? inv(width) : wide/sinh(wide*width)
            pref=BigFloat(ρ)/(2BigFloat(π))
            for id in (:Zhao2020,:Wedepohl1973)
                f=II.Formula(id)(a,b,ρ,mr,s)
                expected=if id===:Zhao2020
                    (pref*mc/a,pref*mc/b,pref*ms/sqrt(BigFloat(a)*b))
                else
                    correction=pref/(BigFloat(b)*(BigFloat(a)+b))
                    (pref*mc/a-correction,pref*mc/b+correction,
                        2pref*ms/(BigFloat(a)+b))
                end
                for (kind,reference) in zip((:inner,:outer,:mutual),expected)
                    actual=f(Val(kind))
                    @test actual isa Complex{T}
                    @test isfinite(actual)
                    @test actual ≈ Complex{T}(reference) rtol=tol atol=eps(T)^2
                end
            end
        end
        rho=T[Inf,100]; epsilon=T[ε,10ε]; mu=T[μ,μ]
        for frequency in T.((1,50,1e5))
            s=complex(zero(T),2T(π)*frequency)
            g=sqrt(s*μ/rho[2]); p=inv(g)
            for (i,k,heights,x) in ((1,1,(T(5),T(5)),T(.02)),
                    (1,2,(T(5),T(8)),T(7)))
                pair=E.EarthPair(i,k,heights,x,(1,1))
                H=sum(heights); lateral=i==k ? zero(T) : x
                direct=i==k ? x : hypot(x,heights[1]-heights[2])
                q=lateral/H; u=H/(2p)
                J=log(((1+inv(u))^2+q^2)/(1+q^2))/2-
                    ((1+u*(1+im*q))^(-3)+(1+u*(1-im*q))^(-3))/24
                reference=s*μ/(2T(π))*(log(hypot(H,lateral)/direct)+J)
                f=EI.Formula(:Alvarado1983)(rho,epsilon,mu,s,nothing)
                @test f(i==k ? Val(:self) : Val(:mutual),pair) ≈ reference rtol=tol
                referenceP=log(hypot(H,lateral)/direct)/(2T(π)*ε)
                classical=EA.Formula(:Ametani2021)(rho,epsilon,mu,s,nothing)
                @test classical(Val(:mutual),pair) ≈ referenceP rtol=tol
            end
            pair=E.EarthPair(1,1,(-T(2),-T(2)),T(.01),(2,2))
            bridges=EI.Formula(:Bridges1995)(rho,epsilon,mu,s,nothing)
            expected=-s*μ/(2T(π))*log(T(17811)/10000*g*pair.separation/2)
            @test bridges(Val(:self),pair) ≈ expected rtol=tol
            @test_throws ArgumentError bridges(Val(:mutual),pair)
            classical=EA.Formula(:Pollaczek1926)(rho,epsilon,mu,s,nothing)
            @test classical(Val(:self),pair)==zero(s)
            @test classical(Val(:mutual),E.EarthPair(1,2,(T(2),-T(2)),T(1),(1,2)))==zero(s)
            for pair in (pair,E.EarthPair(1,2,(-T(2),-T(3)),T(2),(2,2)))
                zi=EI.Formula(:Theethayi2007)(rho,epsilon,mu,s,nothing)(Val(:mutual),pair)
                yi=EA.Formula(:Theethayi2007)(rho,epsilon,mu,s,nothing)(Val(:mutual),pair)
                @test yi ≈ s*zi/(s*μ*(inv(rho[2])+s*epsilon[2])) rtol=tol
            end
        end
    end
    for owner in (EI,EA,II), alias in keys(owner.ALIASES)
        canonical=owner.Formula(alias)
        for selector in (alias,Val(alias))
            explicit=owner.Formula(selector,owner.routes(canonical),owner.assumptions(canonical))
            @test E.formula_id(explicit)===E.formula_id(canonical)
            @test owner.routes(explicit)==owner.routes(canonical)
            @test owner.assumptions(explicit)==owner.assumptions(canonical)
        end
        route=first(values(owner.routes(canonical)))
        old=LineCableModels.FormulaMethod(Val(alias),route.method,route.arguments...)
        key=first(keys(owner.routes(canonical)))
        mapped=owner.Formula(Val(alias),NamedTuple{(key,)}((old,)))
        @test typeof(getproperty(owner.routes(mapped),key)).parameters[1]===E.formula_id(canonical)
    end
end
