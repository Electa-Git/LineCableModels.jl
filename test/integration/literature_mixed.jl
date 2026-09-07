@testitem "Engine / literature assimilation / mixed-source equivalence" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance
    μ0,ε0=4π*1e-7,8.8541878128e-12
    # For lossless air, choose the outgoing +j boundary value on the
    # negative real axis. A signed -0 from k² subtraction is not a
    # different physical radiation condition.
    decayroot(z)=sqrt(complex(real(z),iszero(imag(z)) ? 0.0 : imag(z)))
    for frequency in (1.0,50.0,10000.0), mur in (1.0,3.0),
        Γp in (0.0im,complex(0.002,0.003))
        s=complex(0.0,2π*frequency)
        rho=[Inf,100.0]; ε=ε0.*[1.0,10.0]; μ=μ0.*[1.0,mur]
        κ=inv.(rho).+s.*ε
        # Pawlik's e^(-Γp z) uses k_i^2=-s μ_i κ_i.
        # Engine's longitudinal spectral wave number is k_x=j Γp.
        k2=-s.*μ.*κ
        radial2=Γp^2 .+ k2
        Γengine=im*Γp
        pair=E.EarthPair(1,2,(5.0,-1.5),1.0,(1,2))
        magnetic=quadgk(0.0,Inf;rtol=1e-10) do λ
            u=decayroot.(λ^2 .- radial2)
            exp(-5u[1]-1.5u[2])/(u[1]+μ[1]/μ[2]*u[2])*cos(λ)
        end
        electric=quadgk(0.0,Inf;rtol=1e-10) do λ
            u=decayroot.(λ^2 .- radial2)
            (u[1]+μ[2]/μ[1]*u[2])*exp(-5u[1]-1.5u[2])/
                ((u[1]+μ[1]/μ[2]*u[2])*
                 (k2[2]/k2[1]*u[1]+μ[2]/μ[1]*u[2]))*cos(λ)
        end
        Zpawlik=s*μ[1]/π*first(magnetic)
        Ypawlik=π*κ[1]/first(electric)
        zi=EI.Formula(:MartinsBritto2024)(rho,ε,μ,s,Γengine)
        pa=EA.Formula(:MartinsBritto2024)(rho,ε,μ,s,Γengine)
        @test zi(Val(:mutual),pair) ≈ Zpawlik rtol=3e-8
        # This is a normalization of the single pair, not a matrix inverse.
        @test pa(Val(:mutual),pair) ≈ s/Ypawlik rtol=3e-8
        reversed=E.EarthPair(2,1,(-1.5,5.0),1.0,(2,1))
        @test zi(Val(:mutual),reversed) ≈ Zpawlik rtol=3e-8
        @test pa(Val(:mutual),reversed) ≈ s/Ypawlik rtol=3e-8
        if iszero(Γp)
            # Dawalibi (8), taking γ_i²=s μ_i κ_i so transverse
            # attenuation has inverse-length units.
            gamma2=s.*μ.*κ
            literal=quadgk(0.0,Inf;rtol=1e-10) do λ
                a=sqrt.(λ^2 .+ gamma2)
                exp(-5a[1])*exp(-1.5a[2])/
                    (a[1]*μ[2]+a[2]*μ[1])*cos(λ)
            end
            @test zi(Val(:mutual),pair) ≈
                s*μ[1]*μ[2]/π*first(literal) rtol=3e-8
        end
    end
    for T in (Float32,Float64,BigFloat),layerpair in ((1,1),(2,2),(1,2))
        s=complex(zero(T),100T(π))
        rho=T[Inf,100];epsilon=T(ε0)*T[1,10];mu=T(μ0)*T[1,3]
        heights=layerpair==(1,1) ? (T(5),T(2)) :
            layerpair==(2,2) ? (T(-.5),T(-1.5)) : (T(5),T(-1.5))
        pair=E.EarthPair(1,2,heights,T(1),layerpair)
        value=EA.Formula(:MartinsBritto2024)(
            rho,epsilon,mu,s,complex(zero(T)))(Val(:mutual),pair)
        @test value isa Complex{T}
        @test isfinite(value)
    end
end
