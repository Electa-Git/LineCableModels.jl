@testitem "Engine / literature assimilation / Lucca secondary equation" begin
    using LineCableModels
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0,ε0=4π*1e-7,8.8541878128e-12
    for frequency in (1.0,50.0,10000.0,1e6), resistivity in (10.0,1000.0)
        s=complex(0.0,2π*frequency); γ=sqrt(s*μ0/resistivity)
        f=EI.Formula(:Lucca1994)([Inf,resistivity],[ε0,10ε0],[μ0,μ0],s,nothing)
        for x in (0.0,1.0,30.0), (y1,y2) in ((10.0,-1.0),(2.0,-4.0))
            Y=y1-y2+2/γ; R=sqrt(x^2+Y^2); D=hypot(x,y1-y2)
            reference=s*μ0/(2π)*(log(R/D)-2Y/(3γ^3)*(Y^2-3x^2)/R^6)
            pair=E.EarthPair(1,2,(y1,y2),x,(1,2))
            reversed=E.EarthPair(2,1,(y2,y1),x,(2,1))
            @test f(Val(:mutual),pair) ≈ reference rtol=1e-12
            @test f(Val(:mutual),reversed) ≈ reference rtol=1e-12
        end
    end
end
