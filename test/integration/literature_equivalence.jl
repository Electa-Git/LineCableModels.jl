@testitem "Engine / literature assimilation / equivalent formula witnesses" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance; II=E.InternalImpedance
    μ, ε = 4π*1e-7, 8.8541878128e-12
    rho=[Inf,100.0]; epsilon=[ε,10ε]; permeability=[μ,μ]
    for frequency in (1.0,50.0,1e4)
        s=complex(0.0,2π*frequency)
        for radius in (0.001,0.01), resistivity in (1.724e-8,1e-6)
            # Wait (3b): k_w^2=-s μ σ, with j k_w equal to the
            # positive-real diffusion root. This checks only the good-conductor
            # reduction, not Wait's propagation-dependent parent (3a).
            kw=sqrt(-s*μ/resistivity)
            wait=sqrt(s*μ*resistivity)/(2π*radius) *
                 besselix(0,im*kw*radius)/besselix(1,im*kw*radius)
            schel=II.Formula(:Schelkunoff1934)(0.0,radius,resistivity,1.0,s)(Val(:outer))
            @test schel ≈ wait rtol=1e-11
        end
        for (hi,hj,x) in ((10.0,15.0,3.0),(5.0,5.0,0.2))
            pair=E.EarthPair(1,2,(hi,hj),x,(1,1))
            # Ametani (12)-(13) and Xue (20) use the same difference:
            # m1^2-m0^2 = k_a^2-k_e^2.
            m0sq=s*μ*s*ε; m1sq=s*μ*(0.01+s*10ε)
            kasq=-m0sq; kesq=-m1sq
            @test m1sq-m0sq == kasq-kesq
            integral=quadgk(t->exp(-(hi+hj)*t)*cos(x*t)/
                (t+sqrt(t^2+kasq-kesq)),0.0,Inf;rtol=1e-10)[1]
            witness=s*μ/(2π)*(log(hypot(x,hi+hj)/hypot(x,hi-hj))+2integral)
            wise=EI.Formula(:Wise1934)(rho,epsilon,permeability,s,nothing)(Val(:mutual),pair)
            @test wise ≈ witness rtol=3e-8
        end
        # Uribe (2e), after rationalizing the inherited Pollaczek denominator.
        h_air,h_earth,x=5.0,2.0,1.0
        p=inv(sqrt(s*μ/100)); ξ=h_earth/abs(p); η=x/h_earth; ζ=h_air/h_earth
        integral=quadgk(0.0,Inf;rtol=1e-10) do u
            F=sqrt((u^2+hypot(u^2,1))/2)
            G=inv(2F)
            # F-u=G^2/(F+u) avoids cancellation in Uribe's printed expression.
            (G^2/(F+u)+im*G)*exp(-ξ*(ζ*u+F))*exp(-im*ξ*G)*cos(ξ*η*u)
        end
        uribe=imag(s)*μ/π*integral[1]
        mixed=E.EarthPair(1,2,(h_air,-h_earth),x,(1,2))
        pollaczek=EI.Formula(:Pollaczek1926)(rho,epsilon,permeability,s,nothing)
        @test pollaczek(Val(:mutual),mixed) ≈ uribe rtol=3e-8

        # Nguyen (38)-(41): independent dimensionless normalization.
        pair=E.EarthPair(1,2,(-1.0,-2.0),0.75,(2,2))
        α=imag(s)*μ/100; hm=1.5; x=pair.separation
        J=quadgk(0.0,Inf;rtol=1e-10) do u
            v=sqrt(complex(u^2,1.0))
            hm*sqrt(α)*im/(v+u)*exp(-2hm*sqrt(α)*v)*cos(x*sqrt(α)*u)
        end
        d=hypot(x,1.0); D=hypot(x,3.0)
        nguyen=s*μ/(2π)*(besselk(0,d*sqrt(im*α))-besselk(0,D*sqrt(im*α)))+
            sqrt(imag(s)*μ*100)/(π*hm)*J[1]
        @test pollaczek(Val(:mutual),pair) ≈ nguyen rtol=3e-8
    end
end
