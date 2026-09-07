@testitem "Engine / literature assimilation / Wise permeable-earth reduction" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine
    EI=E.EarthImpedance
    for T in (Float32,Float64,BigFloat), mur in (1,5)
        pi_t=one(T)*π
        μ=T(4)*pi_t/T(10)^7
        ε=T(88541878128)/T(10)^22
        rho=T[Inf,100]; epsilon=T[ε,10ε]; permeability=T[μ,mur*μ]
        s=complex(zero(T),T(100)*pi_t)
        pair=E.EarthPair(1,2,(T(10),T(15)),T(3),(1,1))
        f=EI.Formula(:Wise1931)(rho,epsilon,permeability,s,nothing)
        z=f(Val(:mutual),pair)
        gamma=sqrt(s*permeability[2]/rho[2])
        integral=quadgk(lambda->permeability[2]*exp(-25lambda)*cos(3lambda)/
            (lambda*permeability[2]+sqrt(lambda^2+gamma^2)*μ),
            zero(T),T(Inf);rtol=max(T(1e-9),eps(T)))[1]
        reference=s*μ/(2pi_t)*(log(hypot(T(3),T(25))/hypot(T(3),T(5)))+2integral)
        @test z isa Complex{T}
        @test z ≈ reference rtol=max(T(3e-8),64eps(T))
        @test iszero(f.state.gamma[1])
        @test f.state.gamma[2]^2 ≈ s*permeability[2]/rho[2]
        if mur==1
            carson=EI.Formula(:Carson1926)(rho,epsilon,permeability,s,nothing)(Val(:mutual),pair)
            @test z ≈ carson rtol=max(T(3e-8),64eps(T))
        end
    end
end
