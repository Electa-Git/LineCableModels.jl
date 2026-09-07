@testitem "Engine / literature assimilation / finite-pipe low-frequency limits" begin
    using LineCableModels, SpecialFunctions, LinearAlgebra
    E=LineCableModels.Engine; P=E.PipeImpedance
    μ0=4π*1e-7
    function printed_inner(a,b,mu)
        delta=(b-a)*(b+a)
        return mu/(2π)*(b^4/delta^2*log(b/a)-(3b*b-a*a)/(4delta))
    end
    function parent_harmonic(n,x1,x2,mur)
        ai1=n*(mur+1)/x1*besseli(n,x1)-besseli(n-1,x1)
        ai2=n*(mur-1)/x2*besseli(n,x2)+besseli(n-1,x2)
        ak1=n*(mur+1)/x1*besselk(n,x1)+besselk(n-1,x1)
        ak2=n*(mur-1)/x2*besselk(n,x2)-besselk(n-1,x2)
        return 2mur/x1*(ak2*besseli(n,x1)-ai2*besselk(n,x1))/
            (ai1*ak2-ai2*ak1)
    end
    for T in (Float32,Float64,BigFloat), ratio in (1.000001,1.01,1.3,3.0),
            murvalue in (1.0,3.0,100.0)
        a=T(1)/10; b=a*T(ratio); rho=T(1)/10^6; mur=T(murvalue)
        s=complex(zero(T),T(2)*T(π)/1000)
        f=P.Formula(:Hoidalen2013)(a,b,rho,mur,s)
        dc=rho/(T(π)*(b-a)*(b+a))
        reference=setprecision(BigFloat,256) do
            printed_inner(BigFloat(a),BigFloat(b),4big(π)*big"1e-7"*BigFloat(mur))
        end
        tol=T===Float32 ? 2e-5 : 5e-11
        @test imag(f(Val(:inner)))/imag(s) ≈ reference rtol=tol
        @test real(f(Val(:inner))) == dc
        @test real(f(Val(:outer))) == dc
        @test real(f(Val(:mutual))) == dc
        combination=imag(f(Val(:inner))+f(Val(:outer))-2f(Val(:mutual)))/imag(s)
        @test combination ≈ T(μ0)*mur/(2T(π))*log1p((b-a)/a) rtol=tol
        @test f(Val(:inner)) isa Complex{T}
        pair=P.Pair(1,1,((a/3,zero(T)),(a/3,zero(T))),(a/10,a/10))
        @test f(Val(:self),pair) isa Complex{T}
        @test isfinite(f(Val(:self),pair))
    end
    for mur in (1.0,3.0,100.0), b in (.101,.13,.3), n in (1,2,8)
        a=.1; s=2π*1e-8im; m=sqrt(s*μ0*mur/1e-6)
        @test P._static_pipe_harmonic(n,log(b/a),mur) ≈
            parent_harmonic(n,m*a,m*b,mur) rtol=2e-8
    end
    # Independent full surface Bessel equations (8), (20), (21), with enough
    # precision to retain their tiny reactive parts above the DC resistance.
    setprecision(BigFloat,256) do
        I(n,z)=sum((z/2)^(2k+n)/(factorial(big(k))*factorial(big(k+n))) for k in 0:6)
        function K(n,z)
            @assert abs(z)<big"1e-5"
            # Convergent integer-order K series, independently of production
            # quadrature. Seven terms suffice at these tiny skin arguments.
            harmonic(k)=sum(inv(BigFloat(j)) for j in 1:k)
            q=z^2/4; logarithm=log(z/2)+BigFloat(Base.MathConstants.eulergamma)
            if n==0
                return -logarithm*I(0,z)+
                    sum(harmonic(k)*q^k/factorial(big(k))^2 for k in 1:6)
            end
            return I(0,z)/z+logarithm*I(1,z)-
                2/z*sum(k*harmonic(k)*q^k/factorial(big(k))^2 for k in 1:6)
        end
        for ratio in (big"1.01",big"1.3",big"3.0"), mur in (big"1.0",big"10.0")
            a=big".1"; b=a*ratio; rho=big"1e-6"; mu=4big(π)*big"1e-7"*mur
            s=complex(big"0",2big(π)*big"1e-12"); m=sqrt(s*mu/rho)
            x1,x2=m*a,m*b
            denominator=I(1,x2)*K(1,x1)-I(1,x1)*K(1,x2)
            inner=s*mu/(2big(π)*x1)*(I(0,x1)*K(1,x2)+I(1,x2)*K(0,x1))/denominator
            outer=s*mu/(2big(π)*x2)*(I(0,x2)*K(1,x1)+I(1,x1)*K(0,x2))/denominator
            transfer=s*mu/(2big(π)*x1*x2)/denominator
            f=P.Formula(:Hoidalen2013)(a,b,rho,mur,s)
            for (route,value) in ((Val(:inner),inner),(Val(:outer),outer),(Val(:mutual),transfer))
                @test real(f(route)) ≈ real(value) rtol=big"1e-20"
                @test imag(f(route)) ≈ imag(value) rtol=big"1e-20"
            end
        end
    end
    for mur in (1.0,10.0), frequency in (1e-6,1e-4)
        a=.127; b=.133; s=2π*frequency*im
        f=P.Formula(:Hoidalen2013)(a,b,1e-6,mur,s)
        z=[.04+0im,-.02+.03im,-.025-.02im]; r=[.01,.012,.009]
        pairs=[P.Pair(i,j,((real(z[i]),imag(z[i])),(real(z[j]),imag(z[j]))),(r[i],r[j]))
            for i in 1:3,j in 1:3]
        matrix=[f(i==j ? Val(:self) : Val(:mutual),pairs[i,j]) for i in 1:3,j in 1:3]
        @test matrix ≈ transpose(matrix) rtol=1e-12
        @test minimum(eigvals(Symmetric(real.(matrix)))) >= -1e-15
        @test minimum(eigvals(Symmetric(imag.(matrix)/imag(s)))) > 0
        negative=P.Formula(:Hoidalen2013)(a,b,1e-6,mur,-s)
        @test negative(Val(:mutual),pairs[1,2]) ≈ conj(matrix[1,2]) rtol=1e-12
        if isone(mur)
            for i in 1:3,j in 1:3
                distance=i==j ? r[i] : abs(z[i]-z[j])
                @test matrix[i,j] ≈ f(Val(:inner))+s*μ0/(2π)*log(a/distance) rtol=1e-12
            end
            # Eq. (13) cancels the image logarithm for both self and mutual.
            w=z[1]*conj(z[2])/a^2
            @test sum(P._static_pipe_harmonic(n,log(b/a),mur)*real(w^n) for n in 1:60) ≈
                -log(abs(1-w)) rtol=1e-12
        end
    end
    @test_throws DomainError P.Formula(:Hoidalen2013)(.1,.1,1e-6,1.0,1im*1.0)
    @test_throws DomainError P.Formula(:Hoidalen2013)(.1,.13,1e-6,1.0,0im*1.0)
end
