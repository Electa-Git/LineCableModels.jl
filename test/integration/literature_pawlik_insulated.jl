@testitem "Engine / literature assimilation / Pawlik insulated-wire external terms" begin
    using LineCableModels, SpecialFunctions, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    # Positive-frequency outgoing boundary: u=j*gamma at lambda=0,
    # approached from a medium with positive infinitesimal conductivity.
    outgoing(z)=iszero(imag(z)) && real(z)<0 ? complex(0.,sqrt(-real(z))) : sqrt(z)
    for frequency in (1.,50.,1e5), height in (.6,5.), mur in (1.,5.), source in (1,2)
        s=2π*frequency*im
        rho=[Inf,100.]; epsilon=[ε0,10ε0]; mu=[μ0,mur*μ0]
        k0=imag(s)*sqrt(μ0*ε0)
        for kx in (0.0im,complex(1.2k0+.001,.0001))
            # Source attenuation convention e^-Gamma*z, engine Fourier kx=jGamma.
            Gamma=-im*kx
            k_squared=-s .* mu .* (inv.(rho).+s.*epsilon)
            gamma_squared=k_squared .+ Gamma^2
            gamma=sqrt.(gamma_squared)
            gamma=map(g->imag(g)>0 ? -g : g,gamma)
            b=.0042; B=hypot(2height,b); other=3-source
            Lambda=besselk(0,im*gamma[source]*b)-besselk(0,im*gamma[source]*B)
            roots(lambda)=(outgoing(lambda^2-gamma_squared[source]),
                outgoing(lambda^2-gamma_squared[other]))
            iz=quadgk(0.0,Inf;rtol=1e-10) do t
                lambda=t/(2height); u1,u2=roots(lambda)
                exp(-2height*u1)*cos(lambda*b)/(u1+mu[source]/mu[other]*u2)/(2height)
            end[1]
            iy=quadgk(0.0,Inf;rtol=1e-10) do t
                lambda=t/(2height); u1,u2=roots(lambda)
                numerator=(u1+mu[other]/mu[source]*u2)*exp(-2height*u1)
                denominator=(u1+mu[source]/mu[other]*u2)*
                    (k_squared[other]/k_squared[source]*u1+mu[other]/mu[source]*u2)
                numerator/denominator*cos(lambda*b)/(2height)
            end[1]
            expectedZ=s*mu[source]/(2π)*(Lambda+2iz)
            expectedY=2π*(inv(rho[source])+s*epsilon[source])/(Lambda+2iy)
            h=source==1 ? height : -height
            pair=E.EarthPair(1,1,(h,h),b,(source,source))
            fi=EI.Formula(:Pawlik2018)(rho,epsilon,mu,s,kx)
            fa=EA.Formula(:MartinsBritto2024)(rho,epsilon,mu,s,kx)
            @test fi(Val(:self),pair) ≈ expectedZ rtol=3e-8
            @test fa(Val(:self),pair) ≈ s/expectedY rtol=3e-8
            negative_i=EI.Formula(:Pawlik2018)(rho,epsilon,mu,-s,conj(kx))
            negative_a=EA.Formula(:MartinsBritto2024)(rho,epsilon,mu,-s,conj(kx))
            @test negative_i(Val(:self),pair) ≈ conj(fi(Val(:self),pair)) rtol=1e-11
            @test negative_a(Val(:self),pair) ≈ conj(fa(Val(:self),pair)) rtol=1e-11
            for lambda in (0.,.1,1.,10.)
                u1,u2=roots(lambda)
                literal=(u1+mu[other]/mu[source]*u2)/
                    ((u1+mu[source]/mu[other]*u2)*
                     (k_squared[other]/k_squared[source]*u1+mu[other]/mu[source]*u2))
                @test EA._same_medium_potential_kernel(u1,u2,-k_squared[source],
                    -k_squared[other],mu[source],mu[other]) ≈ literal rtol=1e-12
            end
        end
    end
    # Exact zero transverse constant: neither K0(0)-K0(0) nor a perturbed
    # artificial radius is an admissible numerical representation.
    s=2π*50im; rho=[Inf,100.]; epsilon=[ε0,10ε0]; mu=[μ0,μ0]
    base=EI.Formula(:Pawlik2018)(rho,epsilon,mu,s,nothing)
    kx=sqrt(-base.state.gamma_medium_squared[1])
    fi=EI.Formula(:Pawlik2018)(rho,epsilon,mu,s,kx)
    fa=EA.Formula(:MartinsBritto2024)(rho,epsilon,mu,s,kx)
    @test iszero(fi.state.gamma_medium_squared[1]+fi.state.gamma_squared)
    @test iszero(fa.state.gamma_medium_squared[1]+fa.state.gamma_squared)
    h=.6; b=.0042; B=hypot(2h,b)
    pair=E.EarthPair(1,1,(h,h),b,(1,1))
    z=fi(Val(:self),pair); p=fa(Val(:self),pair)
    @test isfinite(z)
    @test isfinite(p)
    beta_squared=base.state.gamma_medium_squared[2]-base.state.gamma_medium_squared[1]
    ratio=base.state.gamma_medium_squared[2]/base.state.gamma_medium_squared[1]
    iz=quadgk(lambda->exp(-2h*lambda)*cos(b*lambda)/
        (lambda+sqrt(lambda^2+beta_squared)),0.,Inf;rtol=1e-11)[1]
    ip=quadgk(lambda->exp(-2h*lambda)*cos(b*lambda)/
        (ratio*lambda+sqrt(lambda^2+beta_squared)),0.,Inf;rtol=1e-11)[1]
    @test z ≈ s*μ0/(2π)*(log(B/b)+2iz) rtol=3e-8
    @test p ≈ (log(B/b)+2ip)/(2π*ε0) rtol=3e-8
    # Factorization retains the small potential response when the two
    # terms of the older magnetic-plus-coupling expression cancel.
    for N in (1e4,1e8,1e12,1e16)
        a1=complex(1.0); a2=complex(2.0); g1=complex(1.0); g2=complex(N)
        value=EA._same_medium_potential_kernel(a1,a2,g1,g2,1.,1.)
        @test value ≈ inv(N+2) rtol=1e-15
        @test !iszero(value)
    end
end
