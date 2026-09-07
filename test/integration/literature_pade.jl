@testitem "Engine / literature assimilation / De Conti Pade cancellation" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance
    tag=Val(:DeConti2024); μ0=4π*1e-7; ε0=8.8541878128e-12
    function literal(γ,d,H,r)
        D=hypot(H,r); z=γ*D; root=sqrt(1-z)
        i1=(H-8/γ)*r/D^2
        i2=16*(2-z)/z^2*atan(r/(H+D))
        i3=-4*(8-8z+z^2)/(z^2*root)*atan(r*root/(H+D))
        residual=(i1+i2+i3)*exp(-z)
        return E.special_besselk(0,γ*d)+(H^2-r^2)/D^2*E.special_besselk(2,z)-
            2*(H^2-r^2)/(γ^2*D^4)*exp(-γ*H)*(1+γ*H)-2r*H/D^2*residual
    end
    for magnitude in (1e-8,1e-4,0.1,0.49,0.51,0.99,1.01,3.0,30.0),
        angle in (0.4,π/4,1.5), ratio in (0.1,1.0,3.0)
        H=2.0; r=ratio*H; D=hypot(H,r); z=magnitude*cis(angle)
        θmax=atan(r,H)
        expected=exp(-z)*quadgk(θ->-cos(2θ)*(2-z*(cos(θ)-1))/
            (2+z*(cos(θ)-1)),0.0,θmax;rtol=1e-12)[1]
        actual=EI.closed_form_term(tag,Val(:pade_residual),z,H,r,D)
        @test actual ≈ expected rtol=2e-10 atol=1e-14
    end
    for magnitude in (1e-8,0.01,0.49)
        z=BigFloat(magnitude)*cis(BigFloat(π)/4)
        expected=E.special_besselk(2,z)-2/z^2
        actual=EI._sunde_bessel_remainder(ComplexF64(z))
        @test actual ≈ ComplexF64(expected) rtol=1e-12
        for h in (0.1,0.7,1.0)
            expected=(1-(1+h*z)*exp(-h*z))/z^2
            actual=EI._sunde_exponential_remainder(ComplexF64(z),h)
            @test actual ≈ ComplexF64(expected) rtol=1e-12
        end
        γ=z/3; d=BigFloat(0.01); H=BigFloat(2); r=sqrt(BigFloat(5))
        reference=literal(γ,d,H,r)
        actual=EI.closed_form_term(tag,ComplexF64(γ),Float64(d),Float64(H),Float64(r))
        @test actual ≈ ComplexF64(reference) rtol=1e-11
    end
    for frequency in (1.0,50.0,1e4,1e6), resistivity in (10.0,1000.0)
        s=complex(0.0,2π*frequency)
        f=EI.Formula(:DeConti2024)([Inf,resistivity],[ε0,10ε0],[μ0,μ0],s,nothing)
        γ=sqrt(s*μ0*(inv(resistivity)+s*10ε0))
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(-1.0,-1.0),0.02,(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-1.0,-2.0),0.5,(2,2))),
            (Val(:mutual),E.EarthPair(1,2,(-1.0,-2.0),0.0,(2,2))))
            hi,hj=abs.(pair.heights); r=pair.separation; H=hi+hj
            d=hypot(r,hi-hj)
            expected=s*μ0/(2π)*literal(γ,d,H,r)
            @test f(kind,pair) ≈ expected rtol=2e-8
            @test isfinite(f(kind,pair))
        end
    end
    for T in (Float32,Float64,BigFloat)
        s=complex(zero(T),T(100)*T(π))
        f=EI.Formula(:DeConti2024)(T[Inf,100],T[ε0,10ε0],T[μ0,μ0],s,nothing)
        pair=E.EarthPair(1,2,(T(-1),T(-2)),T(0.5),(2,2))
        @test f(Val(:mutual),pair) isa Complex{T}
    end
end
