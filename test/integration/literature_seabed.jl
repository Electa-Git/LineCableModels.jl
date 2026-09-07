@testitem "Engine / literature assimilation / seabed impedance equivalence" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0,ε0=4π*1e-7,8.8541878128e-12
    for frequency in (1.0,50.0,10000.0), hs in (0.1,1.0,10.0),
        relative_mu in ((1.0,1.0,1.0),(1.0,2.0,3.0))
        s=complex(0.0,2π*frequency)
        rho=[Inf,0.2,100.0]; ε=ε0.*[1.0,80.0,10.0]; μ=μ0.*collect(relative_mu)
        gamma=sqrt.(s.*μ.*(inv.(rho).+s.*ε))
        f=EI.Formula(:Tsiamitros2008)(rho,ε,μ,s,nothing,nothing,[Inf,hs,Inf])
        for (di,dj,x,kind) in ((0.5,1.0,0.3,Val(:mutual)),(0.5,0.5,0.02,Val(:self)))
            pair=E.EarthPair(1,kind===Val(:self) ? 1 : 2,(-hs-di,-hs-dj),x,(3,3))
            # Di Lorenzo (30), with h_i=hs+di and h_j=hs+dj.
            correction=quadgk(0.0,Inf;rtol=1e-10) do λ
                a=sqrt.(λ^2 .+ gamma.^2)
                s10=μ[1]*a[2]+μ[2]*a[1]
                d10=μ[1]*a[2]-μ[2]*a[1]
                s21=μ[3]*a[2]+μ[2]*a[3]
                d21=μ[3]*a[2]-μ[2]*a[3]
                e=exp(-2a[2]*hs)
                -(s10*d21-d10*s21*e)/(s10*s21-d10*d21*e)*
                    exp(-a[3]*(di+dj))/a[3]*cos(x*λ)
            end
            source=s*μ[3]/(2π)*(besselk(0,gamma[3]*hypot(x,di-dj))+first(correction))
            @test f(kind,pair) ≈ source rtol=3e-7
            if kind===Val(:mutual)
                reverse=E.EarthPair(2,1,(-hs-dj,-hs-di),x,(3,3))
                @test f(kind,reverse) ≈ source rtol=3e-7
            end
        end
    end
end
