@testitem "Engine / literature assimilation / finite-wall pipe" begin
    using LineCableModels, SpecialFunctions
    const E=LineCableModels.Engine
    const PI=E.PipeImpedance
    function literal_coefficient(n,x1,x2,mur)
        i1,i2=besseli(n,x1),besseli(n,x2)
        k1,k2=besselk(n,x1),besselk(n,x2)
        di1=besseli(n-1,x1)-n/x1*i1
        di2=besseli(n-1,x2)-n/x2*i2
        dk1=-besselk(n-1,x1)-n/x1*k1
        dk2=-besselk(n-1,x2)-n/x2*k2
        ai1=n*mur/x1*i1-di1; ai2=n*mur/x2*i2+di2
        ak1=n*mur/x1*k1-dk1; ak2=n*mur/x2*k2+dk2
        return 2mur/x1*(ak2*i1-ai2*k1)/(ai1*ak2-ai2*ak1)
    end
    @testset "finite-pipe independent source coefficients" begin
        for mur in (1.0,10.0,100.0), frequency in (1.0,50.0,10000.0)
            μ0=4π*1e-7; rho=1e-6; s=complex(0.0,2π*frequency)
            R1,R2=0.127,0.133; m=sqrt(s*μ0*mur/rho)
            f=PI.Formula(:DaSilva2006)(R1,R2,rho,mur,s)
            x1,x2=m*R1,m*R2
            for n in 1:12
                reference=literal_coefficient(n,x1,x2,mur)
                @test PI._finite_harmonic(n,x1,x2,mur) ≈ reference rtol=1e-10
            end
            zpi=s*μ0/(2π)*mur/x1*
                (besseli(0,x1)*besselk(1,x2)+besseli(1,x2)*besselk(0,x1))/
                (besseli(1,x2)*besselk(1,x1)-besseli(1,x1)*besselk(1,x2))
            @test f.state.inner ≈ zpi rtol=1e-11
            for (kind,pair) in (
                (Val(:self),PI.Pair(1,1,((0.05,0.0),(0.05,0.0)),(0.01,0.01))),
                (Val(:mutual),PI.Pair(1,2,((0.05,0.0),(-0.02,0.03)),(0.01,0.01))))
                g=PI._geometry(pair,R1)
                expected=zpi+s*μ0/(2π)*(g.geometric+
                    sum(literal_coefficient(n,x1,x2,mur)*real(g.w^n) for n in 1:40))
                @test f(kind,pair) ≈ expected rtol=1e-9
            end
            center=PI.Pair(1,1,((0.0,0.0),(0.0,0.0)),(0.01,0.01))
            @test f(Val(:self),center) ≈ zpi+s*μ0/(2π)*log(R1/0.01) rtol=1e-11
            reverse=PI.Pair(2,1,((-0.02,0.03),(0.05,0.0)),(0.01,0.01))
            forward=PI.Pair(1,2,((0.05,0.0),(-0.02,0.03)),(0.01,0.01))
            @test f(Val(:mutual),reverse) ≈ f(Val(:mutual),forward) rtol=1e-12
            neg=PI.Formula(:DaSilva2006)(R1,R2,rho,mur,-s)
            @test neg(Val(:mutual),forward) ≈ conj(f(Val(:mutual),forward)) rtol=1e-12
        end
        # A high-order case exercises the wider-range path without truncating the
        # boundary response. At this thickness its infinite-wall limit is exact
        # to far more than double precision.
        n=200; x1=complex(0.01,0.01); x2=1.25x1; mur=1.0
        finite=PI._finite_harmonic(n,x1,x2,mur)
        kprev=besselkx(0,x1)/besselkx(1,x1)
        for k in 1:n-1
            kprev=inv(kprev+2k/x1)
        end
        infinite=2mur/(n*(1+mur)+x1*kprev)
        @test finite ≈ infinite rtol=1e-11
        @test isfinite(finite)
    end
    @testset "pipe precision and geometry" begin
        for T in (Float32,Float64,BigFloat)
            s=complex(zero(T),T(100)*T(π))
            f=PI.Formula(:DaSilva2006)(T(127)/1000,T(133)/1000,T(1)/10^6,one(T),s)
            pair=PI.Pair(1,1,((T(1)/20,zero(T)),(T(1)/20,zero(T))),(T(1)/100,T(1)/100))
            @test f(Val(:self),pair) isa Complex{T}
            @test isfinite(f(Val(:self),pair))
            @test real(f(Val(:self),pair)) > 0
        end
        s=complex(0.0,100π)
        f=PI.Formula(:DaSilva2006;max_terms=1)(0.127,0.133,1e-6,1.0,s)
        @test_throws ErrorException f(Val(:self),PI.Pair(1,1,((0.05,0.0),(0.05,0.0)),(0.01,0.01)))
        @test_throws DomainError f(Val(:self),PI.Pair(1,1,((0.126,0.0),(0.126,0.0)),(0.01,0.01)))
        @test_throws DomainError PI.Formula(:DaSilva2006;rtol=0)
    end

end
