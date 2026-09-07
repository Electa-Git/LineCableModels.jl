@testitem "Engine / literature assimilation / source distance substitutions" begin
    using LineCableModels
    using SpecialFunctions: besselk, hankelh1
    E = LineCableModels.Engine
    EI = E.EarthImpedance
    μ = 4π * 1e-7
    ε = 8.8541878128e-12
    rho, epsilon, permeability = [Inf,100.0], [ε,10ε], [μ,μ]
    s = complex(0.0, 2π*50)
    for x in (0.0,0.75)
        pair = E.EarthPair(1,2,(-1.0,-2.0),x,(2,2))
        d, H = hypot(x,1.0), 3.0
        m = sqrt(s*μ/100)
        gamma = sqrt(s*μ*(0.01+s*10ε))
        euler = Base.MathConstants.eulergamma
        expected_saad = s*μ/(2π)*(besselk(0,m*d)+2exp(-H*m)/(4+m^2*x^2))
        expected_petrache = s*μ/(2π)*log((1+gamma*d)/(gamma*d))
        expected_wedepohl = s*μ/(2π)*(-log(m*d/2)-euler+0.5-2m*H/3)
        expected_vance = imag(s)*μ/(2π*gamma*d) *
                         hankelh1(0,im*gamma*d)/hankelh1(1,im*gamma*d)
        for (id, reference) in (
            :Saad1996=>expected_saad, :Petrache2005=>expected_petrache,
            :Wedepohl1973=>expected_wedepohl, :Vance1978=>expected_vance
        )
            value = EI.Formula(id)(rho,epsilon,permeability,s,nothing)(Val(:mutual),pair)
            @test isfinite(value)
            @test value ≈ reference rtol=1e-12
        end
    end
    for T in (Float32,Float64,BigFloat)
        ss=complex(zero(T), T(2)*(one(T)*π)*T(1e6))
        pair=E.EarthPair(1,2,(-one(T),-T(2)),one(T),(2,2))
        value=EI.Formula(:Vance1978)(T[Inf,0.1],T[ε,10ε],T[μ,μ],ss,nothing)(Val(:mutual),pair)
        @test value isa Complex{T}
        @test isfinite(value)
    end
end

@testitem "Engine / literature assimilation / identifier compatibility" begin
    using LineCableModels
    E=LineCableModels.Engine
    for (owner,old,id) in (
        (E.InternalImpedance,:AmetaniFuse1992,:Ametani1992),
        (E.InternalImpedance,:WedepohlWilcox1973,:Wedepohl1973),
        (E.EarthImpedance,:WedepohlWilcox1973,:Wedepohl1973),
        (E.EarthImpedance,:Papadopoulos2010,:Papadopoulos2010b),
        (E.EarthAdmittance,:Papadopoulos2010,:Papadopoulos2010b),
        (E.EarthImpedance,:Xue2018,:Xue2018b),
        (E.EarthAdmittance,:Xue2018,:Xue2018b)
    )
        @test formula_id(owner.Formula(old)) === id
        @test typeof(owner.Formula(old)) === typeof(owner.Formula(id))
        @test owner.routes(owner.Formula(old)) == owner.routes(owner.Formula(id))
        @test id in owner.formulas()
        @test old ∉ owner.formulas()
    end
end
