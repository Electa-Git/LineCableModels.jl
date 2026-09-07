@testitem "Engine / literature assimilation / scalar Vance admittance" begin
    using LineCableModels, LinearAlgebra
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance
    for T in (Float32,Float64,BigFloat), frequency in (50,100000)
        πT=one(T)*π; μ=T(4)*πT/T(10)^7; ε=T(88541878128)/T(10)^22
        rho=T[Inf,100]; epsilon=T[ε,10ε]; permeability=T[μ,μ]
        s=complex(zero(T),2πT*T(frequency))
        pair=E.EarthPair(1,1,(-one(T),-one(T)),T(1)/100,(2,2))
        ze=EI.Formula(:Vance1978)(rho,epsilon,permeability,s,nothing)(Val(:self),pair)
        coefficient=EA.Formula(:Vance1978)(rho,epsilon,permeability,s,nothing)(Val(:self),pair)
        gamma_squared=s*μ*(inv(rho[2])+s*epsilon[2])
        @test coefficient isa Complex{T}
        @test s/coefficient ≈ gamma_squared/ze rtol=max(T(1e-12),64eps(T))
        assembled=s*inv(reshape([coefficient],1,1))
        @test assembled[1,1] ≈ gamma_squared/ze rtol=max(T(1e-12),64eps(T))
        # A compatible scalar dependency can be selected without copying its
        # impedance mathematics into the admittance family.
        custom=(functor,p)->EI.earth_impedance(Val(:Theethayi2007),Val(:self),functor,p)
        selected=EA.Formula(:Vance1978;impedance=custom)(rho,epsilon,permeability,s,nothing)
        zcustom=EI.Formula(:Theethayi2007)(rho,epsilon,permeability,s,nothing)(Val(:self),pair)
        @test s/selected(Val(:self),pair) ≈ gamma_squared/zcustom rtol=max(T(1e-12),64eps(T))
        other=E.EarthPair(1,2,(-one(T),-one(T)),T(1)/10,(2,2))
        @test_throws ArgumentError selected(Val(:mutual),other)
        @test_throws ArgumentError selected(Val(:self),other)
    end
end
