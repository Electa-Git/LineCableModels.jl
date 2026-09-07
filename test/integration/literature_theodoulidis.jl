@testitem "Engine / literature assimilation / Theodoulidis exact evaluators" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    evaluations=(:finite_integral,:bessel_product,:single_bessel,:hypergeometric,:recursive)
    original_precision=precision(BigFloat)
    for frequency in (1.0,50.0,10000.0), rho in (10.0,1000.0)
        s=complex(0.0,2π*frequency)
        args=([Inf,rho],[ε0,10ε0],[μ0,μ0],s,nothing)
        parent=EI.Formula(:Pollaczek1926)(args...)
        for pair in (E.EarthPair(1,1,(-1.0,-1.0),0.02,(2,2)),
                E.EarthPair(1,2,(-1.0,-2.0),0.5,(2,2)),
                E.EarthPair(1,2,(-1.0,-2.0),5.0,(2,2)),
                E.EarthPair(1,2,(-1.0,-2.0),0.0,(2,2)))
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            expected=parent(kind,pair)
            for evaluation in evaluations
                f=EI.Formula(:Pollaczek1926;evaluation)(args...)
                @test EI.formula_id(f.state.formula)===:Pollaczek1926
                @test f(kind,pair) ≈ expected rtol=3e-8
            end
        end
    end
    for x in (BigFloat(0),BigFloat("0.2"),BigFloat(1),BigFloat(3))
        H=BigFloat(1); k=Complex{BigFloat}(0.5,0.5)
        expected=quadgk(λ->λ*exp(-H*sqrt(λ^2+k^2))*cos(λ*x),
            BigFloat(0),BigFloat(Inf);rtol=BigFloat("1e-30"))[1]
        for evaluation in evaluations
            @test EI._pollaczek_auxiliary(Val(evaluation),k,H,x) ≈ expected rtol=BigFloat("1e-25")
        end
    end
    for T in (Float32,Float64,BigFloat)
        s=complex(zero(T),T(100)*T(π))
        args=(T[Inf,100],T[ε0,10ε0],T[μ0,μ0],s,nothing)
        pair=E.EarthPair(1,2,(T(-1),T(-2)),T(0.5),(2,2))
        for evaluation in evaluations
            @test EI.Formula(:Pollaczek1926;evaluation)(args...)(Val(:mutual),pair) isa Complex{T}
        end
    end
    @test precision(BigFloat)==original_precision
    @test_throws ArgumentError EI.Formula(:Pollaczek1926;evaluation=:unknown)
    args=([Inf,100.0],[ε0,10ε0],[μ0,μ0],complex(0.0,100π),nothing)
    for pair in (E.EarthPair(1,2,(10.0,15.0),3.0,(1,1)),
            E.EarthPair(1,2,(10.0,-1.0),3.0,(1,2)))
        f=EI.Formula(:Pollaczek1926;evaluation=:hypergeometric)(args...)
        @test f(Val(:mutual),pair)==EI.Formula(:Pollaczek1926)(args...)(Val(:mutual),pair)
    end
end
