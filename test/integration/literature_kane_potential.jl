@testitem "Engine / literature assimilation / Kane pipe potential and field terms" begin
    using LineCableModels
    E=LineCableModels.Engine
    for T in (Float32,Float64,BigFloat)
        μ0=T(4)*T(π)/T(10)^7; ε0=T(88541878128)*T(10)^(-22)
        R=T(127)/1000; a,b=T(10)/1000,T(12)/1000
        ri,rj=T(4)/100,T(3)/100
        for εr in T.((1,2.3,10)), θ in T.((0,0.7,2.5))
            ε=ε0*εr
            f=E.PipeAdmittance.Formula(:Kane1995)(R,ε)
            xi=(ri,zero(T)); xj=(rj*cos(θ),rj*sin(θ))
            self=E.PipeImpedance.Pair(1,1,(xi,xi),(a,a))
            # Use a well-separated angle/radius combination for mutual pairs.
            distance=hypot(xi[1]-xj[1],xi[2]-xj[2])
            Lii=μ0/(2T(π))*log((R^2-ri^2)/(R*a))
            Cii=μ0*ε/Lii
            @test f(Val(:self),self) ≈ inv(Cii)
            @test E.PipeImpedance._geometry(self,R).geometric*μ0/(2T(π)) ≈ Lii
            distance>=a+b || continue
            pair=E.PipeImpedance.Pair(1,2,(xi,xj),(a,b))
            argument=(R^4+(ri*rj)^2-2ri*rj*R^2*cos(θ))/
                (R^2*(ri^2+rj^2-2ri*rj*cos(θ)))
            Lji=μ0/(4T(π))*log(argument)
            @test f(Val(:mutual),pair) ≈ Lji/(μ0*ε)
            @test E.PipeImpedance._geometry(pair,R).geometric*μ0/(2T(π)) ≈ Lji
            reversed=E.PipeImpedance.Pair(2,1,(xj,xi),(b,a))
            @test f(Val(:mutual),reversed) ≈ f(Val(:mutual),pair)
            @test f(Val(:mutual),pair) isa T
            # The conducting-dielectric extension is separately composed.
            s=complex(zero(T),T(100)*T(π)); κ=T(1e-5)+s*ε
            P=s*ε/κ*f(Val(:mutual),pair)
            @test P ≈ s/(2T(π)*κ)*log(sqrt(argument))
        end
    end
    @test_throws DomainError E.PipeAdmittance.Formula(:Kane1995)(0.0,1.0)
    @test_throws DomainError E.PipeAdmittance.Formula(:Kane1995)(1.0,0.0)
    base=E.PipeAdmittance.routes(E.PipeAdmittance.Formula(:Kane1995)).self
    custom=(f,p)->2base(f,p)
    selected=Formulation(pipe_admittance=formula(:Kane1995;self=custom))
    @test selected.methods.pipe_admittance.routes.self === custom
end
