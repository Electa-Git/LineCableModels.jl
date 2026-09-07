@testitem "Engine / literature assimilation / Ametani displacement-current image" begin
    using LineCableModels
    E=LineCableModels.Engine; EI=E.EarthImpedance
    @test EI.formula_id(EI.Formula(:Ametani2014)) === :Pettersson1994
    for T in (Float32,Float64,BigFloat), frequency in (1,50,1000000), er in (1,10,100)
        πT=one(T)*π; μ=T(4)*πT/T(10)^7; ε=T(88541878128)/T(10)^22
        rho=T[Inf,100]; epsilon=T[ε,er*ε]; permeability=T[μ,μ]
        s=complex(zero(T),2πT*T(frequency))
        p=inv(sqrt(s*μ*(inv(rho[2])+s*(epsilon[2]-epsilon[1]))))
        f=EI.Formula(:Ametani2014)(rho,epsilon,permeability,s,nothing)
        pairs=(
            (Val(:self),E.EarthPair(1,1,(T(10),T(10)),T(1)/100,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(T(10),T(15)),T(3),(1,1)))
        )
        for (kind,pair) in pairs
            H=sum(pair.heights); x=pair.separation
            ratio=kind === Val(:self) ? (H+2p)/x :
                sqrt((H+2p)^2+x^2)/hypot(x,pair.heights[1]-pair.heights[2])
            @test f(kind,pair) isa Complex{T}
            @test f(kind,pair) ≈ s*μ/(2πT)*log(ratio) rtol=max(T(1e-12),64eps(T))
            if er==1
                conductive=EI.Formula(:Dubanton1969)(rho,epsilon,permeability,s,nothing)
                @test f(kind,pair) ≈ conductive(kind,pair) rtol=max(T(1e-12),64eps(T))
            end
        end
    end
    μ,ε=4π*1e-7,8.8541878128e-12; s=complex(0.0,2π*50)
    @test_throws DomainError EI.Formula(:Ametani2014)([1e9,100.0],[ε,10ε],[μ,μ],s,nothing)
    f=EI.Formula(:Ametani2014)([Inf,100.0],[ε,10ε],[μ,μ],s,nothing)
    @test_throws ArgumentError f(Val(:mutual),E.EarthPair(1,2,(10.0,-1.0),3.0,(1,2)))
end
