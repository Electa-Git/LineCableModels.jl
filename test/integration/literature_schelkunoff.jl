@testitem "Engine / literature assimilation / Schelkunoff boundary currents" begin
    using LineCableModels,SpecialFunctions,LinearAlgebra
    E=LineCableModels.Engine;II=E.InternalImpedance
    rho=1.7241e-8;mu0=4π*1e-7
    # Solve the original (71) for both independent current excitations.
    function boundary_solution(a,b,rho,mur,s)
        m=sqrt(s*mu0*mur/rho)
        field=[besseli(1,m*a) besselk(1,m*a);
            besseli(1,m*b) besselk(1,m*b)]
        coefficients=field\Diagonal([-1/(2π*a),1/(2π*b)])
        electric=rho*m*[besseli(0,m*a) -besselk(0,m*a);
            besseli(0,m*b) -besselk(0,m*b)]
        return electric*coefficients
    end
    for a in (.001,.008,.0099),mur in (1.,10.),f in (1e-6,50.,1e4)
        s=2π*f*im;b=.01
        functor=II.Formula(:Schelkunoff1934)(a,b,rho,mur,s)
        surfaces=[functor(Val(:inner)) functor(Val(:mutual));
            functor(Val(:mutual)) functor(Val(:outer))]
        @test surfaces ≈ boundary_solution(a,b,rho,mur,s) rtol=2e-10
        @test minimum(eigvals(Symmetric(real.(surfaces))))>=-1e-12
        reverse=II.Formula(:Schelkunoff1934)(a,b,rho,mur,-s)
        for kind in (:inner,:outer,:mutual)
            @test reverse(Val(kind)) ≈ conj(functor(Val(kind))) rtol=1e-12
        end
    end
    for T in (Float32,Float64,BigFloat),a in (0.,.005,1e-8),f in (0.,1e-8,50.)
        b=T(.01);inner=T(a);resistivity=T(rho)
        s=complex(zero(T),T(2)*T(π)*T(f))
        functor=II.Formula(:Schelkunoff1934)(inner,b,resistivity,one(T),s)
        for kind in (:inner,:outer,:mutual)
            value=functor(Val(kind))
            @test value isa Complex{T}
            @test isfinite(value)
            if iszero(f)
                expected=iszero(a)&&kind!==:outer ? zero(T) :
                    resistivity/(T(π)*(b^2-inner^2))
                @test value ≈ expected rtol=10eps(T)
            end
        end
        if f==1e-8
            @test functor(Val(:outer)) ≈ resistivity/(T(π)*(b^2-inner^2)) rtol=2e-6
        end
    end
    # The old intermediate exp(+skin_depths) overflowed in this regime.
    for f in (1e8,1e12),mur in (1.,100.)
        s=2π*f*im;a=.005;b=.01
        functor=II.Formula(:Schelkunoff1934)(a,b,rho,mur,s)
        skin=sqrt(s*mu0*mur*rho)/(2π)
        @test functor(Val(:outer)) ≈ skin/b rtol=.001
        @test functor(Val(:inner)) ≈ skin/a rtol=.001
        @test abs(functor(Val(:mutual)))<=1e-100
        @test isfinite(functor(Val(:outer))) && isfinite(functor(Val(:inner)))
    end
    for args in ((.01,.01,rho,1.),(-.01,.02,rho,1.),(0.,.01,-rho,1.))
        @test_throws DomainError II.Formula(:Schelkunoff1934)(args...,100π*im)
    end
end
