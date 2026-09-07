@testitem "Engine / literature assimilation / Lee multilayer parent" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EI=E.EarthImpedance
    μ0,ε0=4π*1e-7,8.8541878128e-12
    # Original equations (4)-(9), printed p. 1382. F and G use absolute
    # interface depths, not the local layer thickness used by the engine.
    function source_kernel(λ,gamma_squared,air_squared,depths)
        u=sqrt.(λ^2 .+ gamma_squared .- air_squared)
        n=length(u)
        n==1 && return 2/(λ+only(u))
        t=cumsum(depths)
        F=(u[n-1]+u[n])/μ0
        G=(u[n-1]-u[n])*exp(-2u[n-1]*t[n-1])/μ0
        for m in (n-2):-1:1
            Fnext=((u[m]+u[m+1])*F+
                (u[m]-u[m+1])*G*exp(2u[m+1]*t[m]))/μ0
            Gnext=((u[m]-u[m+1])*F+
                (u[m]+u[m+1])*G*exp(2u[m+1]*t[m]))*exp(-2u[m]*t[m])/μ0
            F,G=Fnext,Gnext
        end
        return 2*(F+G)/((λ+u[1])*F+(λ-u[1])*G)
    end
    for count in 1:4, frequency in (1.0,50.0,10000.0,1e6)
        s=complex(0.0,2π*frequency); air_squared=s^2*μ0*ε0
        soilrho=[100.0,500.0,50.0,200.0][1:count]
        soileps=[10.0,6.0,20.0,8.0][1:count]
        depths=[2.0,3.0,5.0][1:count-1]
        rho=[Inf;soilrho]; ε=ε0*[1.0;soileps]; μ=fill(μ0,count+1)
        gamma_squared=[s*μ0*(inv(r)+s*ε0*e) for (r,e) in zip(soilrho,soileps)]
        Γ=sqrt(-air_squared)
        f=EI.Formula(:Lee2014)(rho,ε,μ,s,nothing,nothing,[Inf;depths;Inf])
        for x in (0.0,3.0,20.0)
            pair=E.EarthPair(1,2,(10.0,15.0),x,(1,1))
            # At H=25 m, the omitted lambda>5 tail is below 10^-50.
            correction=quadgk(λ->source_kernel(λ,gamma_squared,air_squared,depths)*
                exp(-25λ)*cos(x*λ),0.0,5.0;rtol=1e-10)[1]
            source=s*μ0/(2π)*(log(hypot(x,25.0)/hypot(x,5.0))+correction)
            @test f(Val(:mutual),pair) ≈ source rtol=3e-8
        end
    end
end
