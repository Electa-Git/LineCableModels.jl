@testitem "Engine / literature assimilation / guarded complex Bessel evaluations" begin
    using LineCableModels,SpecialFunctions,QuadGK
    E=LineCableModels.Engine
    for n in (0,1,12),x in (1e-8,.1,1.,10.,100.,300.),sense in (-1,1)
        z=complex(x,sense*x)
        result=setprecision(BigFloat,128) do
            E.special_besselix(n,Complex{BigFloat}(z))
        end
        @test result ≈ besselix(n,z) rtol=3e-13
    end
    for x in (1e-8,.001,.1,.5),sense in (-1,1),n in (0,1)
        z=complex(x,sense*x)
        result=setprecision(BigFloat,128) do
            E.special_besselk(n,Complex{BigFloat}(z))
        end
        reference=quadgk(t->exp(-z*cosh(t))*cosh(n*t),0.,30.;rtol=1e-12)[1]
        @test result ≈ reference rtol=2e-11
    end
    setprecision(BigFloat,256) do
        for x in (big".000001",big".1",big".5")
            z=x+im*x
            i0=E.special_besselix(0,z)*exp(real(z))
            i1=E.special_besselix(1,z)*exp(real(z))
            k0=E.special_besselk(0,z); k1=E.special_besselk(1,z)
            @test i0*k1+i1*k0 ≈ inv(z) rtol=big"1e-70"
        end
    end
    for n in (0,1,3),z in (-.1+.3im,-10+2im,-30+143im),sense in (-1,1)
        argument=complex(real(z),sense*imag(z))
        result=E.special_besselkx(n,Complex{BigFloat}(argument))
        @test result ≈ besselkx(n,argument) rtol=3e-12
        @test E.special_besselk(n,Complex{BigFloat}(argument))*exp(argument) ≈
            result rtol=3e-12
    end
    @test_throws DomainError E.special_besselkx(1,Complex{BigFloat}(-1,0))
    original=precision(BigFloat)
    # Julia 1.12 precision scopes must not change another task's arithmetic.
    ready=Channel{Nothing}(1); release=Channel{Nothing}(1)
    worker=@async setprecision(BigFloat,original+111) do
        put!(ready,nothing);take!(release)
        @test precision(BigFloat)==original+111
        E.special_besselix(1,Complex{BigFloat}(100+100im))
    end
    take!(ready)
    @test precision(BigFloat)==original
    put!(release,nothing);fetch(worker)
    @test precision(BigFloat)==original
end
