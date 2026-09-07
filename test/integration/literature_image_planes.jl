@testitem "Engine / literature assimilation / complex return planes" begin
    using LineCableModels, LinearAlgebra
    E = LineCableModels.Engine
    EI = E.EarthImpedance
    for T in (Float32, Float64, BigFloat), frequency in (1, 50, 100000)
        πT = one(T) * π
        μ, ε = T(4) * πT / T(10)^7, T(88541878128) / T(10)^22
        s = complex(zero(T), 2πT * T(frequency))
        rho, epsilon, permeability = T[Inf, 100], T[ε, 10ε], T[μ, μ]
        single = E.EarthPair(1, 1, (T(10), T(10)), T(1)/100, (1, 1))
        mutual = E.EarthPair(1, 2, (T(10), T(15)), T(3), (1, 1))
        reversed = E.EarthPair(2, 1, reverse(mutual.heights), mutual.separation, (1, 1))
        p = inv(sqrt(s*μ/100))
        expected_self = s*μ/(2πT)*log(2*(10+p)/single.separation)
        expected_mutual = s*μ/(2πT)*log(sqrt((25+2p)^2+9)/sqrt(T(34)))
        f = EI.Formula(:Dubanton1969)(rho, epsilon, permeability, s, nothing)
        tolerance = max(T(1e-12), 64eps(T))
        for (kind, pair, reference) in (
            (Val(:self), single, expected_self), (Val(:mutual), mutual, expected_mutual)
        )
            @test f(kind, pair) isa Complex{T}
            @test f(kind, pair) ≈ reference rtol=tolerance
            @test EI.Formula(:Deri1981)(rho,epsilon,permeability,s,nothing)(kind,pair) ≈ reference rtol=tolerance
        end
        @test f(Val(:mutual),reversed) ≈ expected_mutual rtol=tolerance
        @test EI.Formula(:Dubanton1969)(rho,epsilon,permeability,-s,nothing)(Val(:self),single) ≈ conj(expected_self) rtol=tolerance
        for old in (:Gary1976,:DeriSemlyen1981)
            @test formula_id(EI.Formula(old)) === :Dubanton1969
            @test old ∉ EI.formulas()
            @test EI.Formula(old)(rho,epsilon,permeability,s,nothing)(Val(:self),single) == f(Val(:self),single)
        end
    end

    μ, ε = 4π*1e-7, 8.8541878128e-12
    pair = E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))
    for frequency in (1.0,50.0,1e4), soils in ([100.0,100.0,100.0],[10.0,100.0,1000.0],[1000.0,100.0,10.0])
        s=complex(0.0,2π*frequency)
        rho=vcat(Inf,soils); epsilon=fill(10ε,4); epsilon[1]=ε; permeability=fill(μ,4)
        thickness=[Inf,0.25,0.5,Inf]
        f=EI.Formula(:Deri1981)(rho,epsilon,permeability,s,nothing,nothing,thickness)
        # Independent downward field-transfer matrices from dE/dx=-s μ H,
        # dH/dx=-σ E, with E/H=ζ in the bottom half-space.
        matrix=Matrix{ComplexF64}(I,2,2)
        for k in 2:3
            g=sqrt(s*μ/rho[k]); ζ=s*μ/g; t=g*thickness[k]
            layer=[cosh(t) -ζ*sinh(t); -sinh(t)/ζ cosh(t)]
            matrix=layer*matrix
        end
        ζbottom=sqrt(s*μ*rho[end])
        E0=(ζbottom*matrix[2,2]-matrix[1,2])/(matrix[1,1]-ζbottom*matrix[2,1])
        p=E0/(s*μ)
        @test f.state.return_plane_depth ≈ p rtol=1e-11
        @test f(Val(:mutual),pair) ≈ s*μ/(2π)*log(sqrt((25+2p)^2+9)/sqrt(34)) rtol=1e-11
        if all(==(first(soils)),soils)
            homogeneous=EI.Formula(:Dubanton1969)(rho[1:2],epsilon[1:2],permeability[1:2],s,nothing)
            @test f(Val(:mutual),pair) ≈ homogeneous(Val(:mutual),pair) rtol=1e-11
        end
    end
    s=complex(0.0,2π*50)
    rho=[Inf,10.0,1000.0]; epsilon=[ε,10ε,10ε]; permeability=fill(μ,3)
    for thickness in ([Inf,0.0,Inf],[Inf,1e9,Inf])
        f=EI.Formula(:Deri1981)(rho,epsilon,permeability,s,nothing,nothing,thickness)
        expected_rho=iszero(thickness[2]) ? rho[3] : rho[2]
        @test f.state.return_plane_depth ≈ inv(sqrt(s*μ/expected_rho)) rtol=1e-12
        @test isfinite(f(Val(:mutual),pair))
    end
    @test_throws DomainError EI.Formula(:Deri1981)(rho,epsilon,permeability,s,nothing,nothing,[Inf,-1.0,Inf])
    @test_throws DimensionMismatch EI.Formula(:Deri1981)(rho,epsilon,permeability,s,nothing)
    @test_throws ArgumentError EI.Formula(:Deri1981)(rho,epsilon,permeability,s,complex(0.01),nothing,[Inf,1.0,Inf])
end

@testitem "Engine / literature assimilation / Saad small-argument variant" begin
    using LineCableModels
    E=LineCableModels.Engine; EI=E.EarthImpedance
    for T in (Float32,Float64,BigFloat)
        πT=one(T)*π; μ=T(4)*πT/T(10)^7; ε=T(88541878128)/T(10)^22
        rho=T[Inf,100]; epsilon=T[ε,10ε]; permeability=T[μ,μ]
        s=complex(zero(T),2πT*T(50)); m=sqrt(s*μ/100)
        f=EI.Formula(:Saad1996;approximation=:small_argument)(rho,epsilon,permeability,s,nothing)
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(-one(T),-one(T)),T(1)/100,(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-one(T),-T(2)),T(3),(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-one(T),-T(2)),zero(T),(2,2)))
        )
            d=hypot(pair.separation,pair.heights[1]-pair.heights[2])
            H=sum(abs,pair.heights)
            expected=s*μ/(2πT)*(-log(m*d/2)-one(T)*Base.MathConstants.eulergamma+T(1)/2-m*H/2)
            @test f(kind,pair) isa Complex{T}
            @test f(kind,pair) ≈ expected rtol=max(T(1e-12),64eps(T))
            parent=EI.Formula(:Saad1996)(rho,epsilon,permeability,s,nothing)(kind,pair)
            @test abs(f(kind,pair)-parent) < T(0.005)*abs(parent)
        end
        @test formula_id(EI.Formula(:Saad1996;approximation=:small_argument)) === :Saad1996
        @test EI.routes(EI.Formula(:Saad1996)) == EI.routes(EI.Formula(:Saad1996;approximation=:closed_form))
    end
    @test_throws ArgumentError EI.Formula(:Saad1996;approximation=:invalid)
    @test_throws ArgumentError EI.Formula(:Saad1996;approximation=:small_argument,unknown=true)
end
