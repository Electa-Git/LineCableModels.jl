@testitem "Engine / literature assimilation / De Conti approximations" begin
    using LineCableModels
    using SpecialFunctions: besselk
    E = LineCableModels.Engine
    EI, EA = E.EarthImpedance, E.EarthAdmittance
    for T in (Float32, Float64, BigFloat)
        pi_t = one(T) * π
        μ = T(4) * pi_t / T(10)^7
        ε = T(88541878128) / T(10)^22
        rho = T[Inf, 100]
        epsilon = T[ε, 10ε]
        permeability = T[μ, μ]
        for frequency in T.((1, 10000, 1000000)), self in (true, false)
            s = complex(zero(T), 2pi_t * frequency)
            pair = self ? E.EarthPair(1, 1, (-one(T), -one(T)), T(1)/100, (2,2)) :
                   E.EarthPair(1, 2, (-one(T), -T(2)), T(3)/4, (2,2))
            a = s * sqrt(μ * ε)
            b = sqrt(s * μ * (inv(rho[2]) + s * epsilon[2]))
            H = sum(abs, pair.heights)
            d = hypot(pair.separation, pair.heights[1] - pair.heights[2])
            D = hypot(pair.separation, H)
            image = (b-a)/(b+a) * exp(-H*b) * 2/(4+b^2*pair.separation^2)
            alpha = (b^2-a^2)/(b^2+a^2)
            euler = one(T) * Base.MathConstants.eulergamma
            direct_small = -log(b*d/2)-euler
            direct_full = E.special_besselk(0, b*d)
            image_full = E.special_besselk(0, b*D)
            for (owner, id, expected) in (
                (EI, :DeConti2023a, s*μ/(2pi_t)*(direct_full+image)),
                (EI, :DeConti2023b, s*μ/(2pi_t)*(direct_small+image)),
                (EA, :DeConti2023a, s/(2pi_t*(inv(rho[2])+s*epsilon[2])) *
                    (direct_full+alpha*image_full)),
                (EA, :DeConti2023b, s/(2pi_t*(inv(rho[2])+s*epsilon[2])) *
                    (log(D/d)-(alpha+1)*(euler+log(b*D/2))))
            )
                functor = owner.Formula(id)(rho, epsilon, permeability, s, nothing)
                z = functor(self ? Val(:self) : Val(:mutual), pair)
                @test z isa Complex{T}
                @test isfinite(z)
                @test z ≈ expected rtol = max(T(2e-12), T(64)*eps(T))
                reversed = E.EarthPair(pair.column, pair.row, reverse(pair.heights),
                    pair.separation, reverse(pair.layers))
                @test functor(Val(:mutual), reversed) ≈ z
                @test_throws ArgumentError functor(Val(:mutual),
                    E.EarthPair(1,2,(one(T),one(T)),one(T),(1,1)))
                @test_throws ArgumentError owner.Formula(id)(
                    rho, epsilon, permeability, s, one(s))
            end
        end
    end

    # The small-argument formulas approach their respective Bessel parents;
    # they are kept as distinct approximations, not deduplicated as identities.
    rho, epsilon, permeability = [Inf,100.0], 8.8541878128e-12 .* [1.0,10.0], fill(4π*1e-7,2)
    pair = E.EarthPair(1,2,(-1.0,-2.0),0.5,(2,2))
    for owner in (EI,EA)
        s=complex(0.0,2π*0.01)
        exact = owner.Formula(:DeConti2023a)(rho,epsilon,permeability,s,nothing)(Val(:mutual),pair)
        small = owner.Formula(:DeConti2023b)(rho,epsilon,permeability,s,nothing)(Val(:mutual),pair)
        @test small ≈ exact rtol=1e-6
    end
end
