@testitem "Engine / literature assimilation / pipe hybrid and finite boundary witnesses" begin
    using LineCableModels, SpecialFunctions, LinearAlgebra
    E=LineCableModels.Engine; P=E.PipeImpedance
    μ0=4π*1e-7
    function fortin_boundary(n,x1,ratio,mur)
        x2=x1*ratio
        di(n,x)=besseli(n-1,x)-n/x*besseli(n,x)
        dk(n,x)=-besselk(n-1,x)-n/x*besselk(n,x)
        # Fourier amplitudes of A and mu^-1 dA/dr are continuous at both
        # radii, source (6)–(7). Radius is scaled by the inner pipe radius.
        matrix=ComplexF64[
            besseli(n,x1) besselk(n,x1) -1 0
            x1/mur*di(n,x1) x1/mur*dk(n,x1) -n 0
            besseli(n,x2) besselk(n,x2) 0 -ratio^(-n)
            x1/mur*di(n,x2) x1/mur*dk(n,x2) 0 n*ratio^(-n-1)
        ]
        scales=vec(maximum(abs.(matrix);dims=1))
        solution=((matrix ./ transpose(scales)) \ ComplexF64[1/n,-1,0,0]) ./ scales
        # The cavity logarithm already includes the perfect-conductor image,
        # so its -1/n harmonic must be restored before adding this response.
        return solution[3]+1/n
    end
    for mur in (1.0,10.0), frequency in (5.0,50.0,500.0),
            thickness_ratio in (1.04,1.5), n in 1:12
        a=.127; b=a*thickness_ratio; rho=1e-6
        m=sqrt(2π*frequency*im*μ0*mur/rho)
        @test P._finite_harmonic(n,m*a,m*b,mur) ≈
            fortin_boundary(n,m*a,thickness_ratio,mur) rtol=1e-10 atol=1e-13
    end
    for mur in (1.0,10.0), frequency in (1.0,50.0,10000.0)
        a=.127; b=.133; rho=1e-6; s=2π*frequency*im
        f=P.Formula(:Yang2001)(a,b,rho,mur,s)
        x=a*sqrt(s*μ0*mur/rho)
        inner=E.InternalImpedance.Formula(:Schelkunoff1934)(a,b,rho,mur,s)(Val(:inner))
        for n in (1,2,8,40)
            literal=2mur/(n*(1+mur)+x*besselk(n-1,x)/besselk(n,x))
            @test P._infinite_pipe_harmonic(n,x,mur) ≈ literal rtol=1e-12
        end
        z=[.04+0im,-.02+.03im]; radii=[.01,.012]
        for i in 1:2,j in 1:2
            pair=P.Pair(i,j,((real(z[i]),imag(z[i])),(real(z[j]),imag(z[j]))),
                (radii[i],radii[j]))
            w=z[i]*conj(z[j])/a^2
            geometric=i==j ? log(a/radii[i])+log(1-abs2(z[i])/a^2) :
                log(abs(a^2-z[i]*conj(z[j]))/(a*abs(z[i]-z[j])))
            # Yang (9)–(12), with dimensionless C_n and one SI prefactor.
            # De Silva (A2) has the identical K_(n-1)/K_n denominator.
            harmonic=sum(2mur*real(w^n)/(n*(1+mur)+
                x*besselk(n-1,x)/besselk(n,x)) for n in 1:60)
            expected=inner+s*μ0/(2π)*(geometric+harmonic)
            @test f(i==j ? Val(:self) : Val(:mutual),pair) ≈ expected rtol=1e-10
        end
        # The thick-wall limit of the complete finite response agrees.
        full=P.Formula(:DaSilva2006)(a,a+40real(inv(sqrt(s*μ0*mur/rho))),rho,mur,s)
        for n in 1:12
            @test P._finite_harmonic(n,full.state.x1,full.state.x2,mur) ≈
                P._infinite_pipe_harmonic(n,x,mur) rtol=1e-10
        end
    end
    for T in (Float32,Float64,BigFloat)
        a=T(127)/1000; b=T(133)/1000; rho=T(1)/10^6
        f=P.Formula(:Yang2001)(a,b,rho,one(T),complex(zero(T),T(100)*T(π)))
        pair=P.Pair(1,1,((a/3,zero(T)),(a/3,zero(T))),(a/10,a/10))
        @test f(Val(:self),pair) isa Complex{T}
        @test isfinite(f(Val(:self),pair))
    end
    # This hybrid must not be identified with a finite magnetic wall at LF.
    a=.127; b=.128; rho=1e-6; mur=100.0; s=2π*1e-4im
    x=a*sqrt(s*μ0*mur/rho)
    @test !isapprox(P._infinite_pipe_harmonic(1,x,mur),
        P._finite_harmonic(1,x,x*b/a,mur);rtol=1e-2)
end
