@testitem "Engine / literature assimilation / Iracheta residual and cutoff" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    for magnitude in (0.1,1.0,20.0,29.9,30.1,50.0), ratio in (0.1,1.0,10.0)
        H=BigFloat(1); x=BigFloat(ratio); R=hypot(H,x)
        k=BigFloat(magnitude)/R*complex(inv(sqrt(BigFloat(2))),inv(sqrt(BigFloat(2))))
        base=(H^2-x^2)/R^4*exp(-k*H)*(1+k*H)
        expected=if magnitude>30
            base
        else
            θ=atan(x,H)
            # Direct source-defined finite moment combination, t=cos(theta).
            moment=quadgk(t->(2sin(t)^2-1)*exp(-k*R*cos(t)),
                BigFloat(0),θ;rtol=BigFloat("1e-30"))[1]
            base+k^2*x*H/R^2*moment
        end
        @test EI._pollaczek_auxiliary(Val(:recursive),k,H,x) ≈ expected rtol=BigFloat("1e-22")
        if magnitude>30
            # The omitted residual gives the normalized Bessel/image formula,
            # not a claim of equality to the full integral in this range.
            ω=Float64(abs(k)^2)*100/μ0; s=complex(0.0,ω)
            pair=E.EarthPair(1,2,(-1/3,-2/3),ratio,(2,2))
            f=EI.Formula(:Pollaczek1926;evaluation=:recursive)(
                [Inf,100.0],[ε0,10ε0],[μ0,μ0],s,nothing)
            g=sqrt(s*μ0/100); distance=hypot(ratio,1/3); image=hypot(ratio,1.0)
            K=g*image
            bracket=besselk(0,g*distance)+(1-ratio^2)/image^2*
                (besselk(0,K)+2besselk(1,K)/K-2*(1+g)*exp(-g)/K^2)
            @test f(Val(:mutual),pair) ≈ s*μ0/(2π)*bracket rtol=1e-10
        end
    end
end
