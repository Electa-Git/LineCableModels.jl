@testitem "Engine / literature assimilation / Noda and Theethayi source substitutions" begin
    using LineCableModels
    E=LineCableModels.Engine;EI=E.EarthImpedance
    mu0=4π*1e-7;eps0=8.8541878128e-12
    for T in (Float32,Float64,BigFloat),freq in (50,100000,10000000),
            soilrho in (10,1000)
        s=complex(zero(T),2T(π)*T(freq))
        rho=T[Inf,soilrho];epsilon=T[eps0,10eps0];mu=T[mu0,mu0]
        mu_source=EI.vacuum_permeability(one(T))
        for id in (:Noda2006,:Theethayi2007)
            formula=EI.Formula(id);leaf=formula(rho,epsilon,mu,s,nothing)
            negative=formula(rho,epsilon,mu,-s,nothing)
            for self in (true,false),x in (T(.02),T(2),T(40),T(100))
                h1,h2=id==:Noda2006 ? (T(10),T(12)) : (T(-.5),T(-1))
                heights=self ? (h1,h1) : (h1,h2)
                layer=id==:Noda2006 ? 1 : 2
                pair=E.EarthPair(1,self ? 1 : 2,heights,x,(layer,layer))
                H=sum(abs,heights)
                if id==:Noda2006
                    lateral=self ? zero(T) : x
                    theta=atan(lateral/H)*180/T(π)
                    A=theta<=T(50.45) ? T(.07360) : T(.002474)*theta-T(.05127)
                    alpha=theta<=T(50.45) ? T(.15) : T(.004726)*theta-T(.08852)
                    beta=(1-A*alpha)/(1-A)
                    p=inv(sqrt(s*mu_source/rho[2]))
                    d=self ? x : hypot(x,h1-h2)
                    expected=s*mu[1]/(2T(π))*
                        (A*log(sqrt((H+2alpha*p)^2+lateral^2)/d)+
                        (1-A)*log(sqrt((H+2beta*p)^2+lateral^2)/d))
                else
                    gamma=sqrt(s*mu_source*(1/rho[2]+s*epsilon[2]))
                    # Source p. 753: horizontal distance and average depth,
                    # even when the two wire depths differ.
                    expected=s*mu[1]/(2T(π))*
                        (log((1+gamma*x)/(gamma*x))+
                        2exp(-H*abs(gamma))/(4+gamma^2*x^2))
                end
                kind=self ? Val(:self) : Val(:mutual)
                value=leaf(kind,pair)
                @test value isa Complex{T}
                @test value ≈ expected rtol=100eps(T)
                @test negative(kind,pair) ≈ conj(value) rtol=100eps(T)
            end
        end
    end
    leaf=EI.Formula(:Theethayi2007)([Inf,100.],[eps0,10eps0],
        [mu0,mu0],100π*im,nothing)
    @test_throws DomainError leaf(Val(:mutual),E.EarthPair(1,2,(-1.,-2.),0.,(2,2)))
    @test_throws ArgumentError leaf(Val(:mutual),E.EarthPair(1,2,(1.,2.),1.,(1,1)))
end
