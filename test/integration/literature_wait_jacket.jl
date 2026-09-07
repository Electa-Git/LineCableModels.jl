@testitem "Engine / literature assimilation / Wait thin jacket modal decomposition" begin
    using LineCableModels
    E=LineCableModels.Engine
    for T in (Float32,Float64,BigFloat), frequency in T.((1,50,1000000)),
        resistivity in T.((100,Inf))
        μ0=T(4)*T(π)/T(10)^7; ε0=T(88541878128)*T(10)^(-22)
        s=complex(zero(T),T(2)*T(π)*frequency)
        ε=ε0*T(3)+(isinf(resistivity) ? zero(s) : inv(resistivity)/s)
        material=Material(:insulator,resistivity,T(3),one(T),T(20),zero(T))
        a,b=T(1)/100,T(12)/1000
        κ=E.InsulationAdmittance.Formula(:Ametani2004)(material,frequency,T(20))
        P=E.potential_coefficient(a,b,κ,s)
        Z=E.InsulationImpedance.Formula(:Ametani1980)(a,b,one(T),s)
        for β in (zero(s),complex(T(0.1),T(-0.02)),complex(T(0.2),T(-0.1)))
            kc2=-s^2*μ0*ε
            source=(β^2-kc2)/(T(2)*T(π)*s*ε)*log(b/a)
            @test Z+β^2*P/s ≈ source rtol=T===Float32 ? T(1e-5) : T(1e-12)
        end
        @test abs(sqrt(-s^2*μ0*ε)*b)<T(0.1)
    end
end
