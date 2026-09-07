@testitem "Engine / literature assimilation / Xue EHEM potential normalization" begin
    using LineCableModels, QuadGK
    E=LineCableModels.Engine; EP=LineCableModels.EarthProps; EH=EP.EHEM
    μ0,ε0=4π*1e-7,8.8541878128e-12
    model=build(EP.EarthModel,(
        EP.EarthLayer(100.0,10.0,1.0,2.0),
        EP.EarthLayer(500.0,6.0,1.0,10.0),
        EP.EarthLayer(50.0,20.0,1.0)))
    rho=collect(getfield.(model.layers,:rho)); eps_r=collect(getfield.(model.layers,:eps_r))
    mu_r=collect(getfield.(model.layers,:mu_r))
    for frequency in (1.0,50.0,1e4,1e6)
        s=complex(0.0,2π*frequency); gamma0=s^2*μ0*ε0
        q=[sqrt(s*μ0*(inv(rho[k])+s*ε0*eps_r[k])-gamma0) for k in 2:4]
        equivalent=q[end]
        for k in 2:-1:1
            g=q[k]; e=exp(-2model.layers[k+1].thickness*g)
            equivalent=g*(g+equivalent-(g-equivalent)*e)/
                (g+equivalent+(g-equivalent)*e)
        end
        for (kind,pair) in (
            (Val(:self),E.EarthPair(1,1,(10.0,10.0),0.02,(1,1))),
            (Val(:mutual),E.EarthPair(1,2,(10.0,15.0),3.0,(1,1))))
            material=EH.Formula(:Xue2021)(
                Val(:overhead),rho,eps_r,mu_r,model,pair,frequency)
            recovered=s*μ0*(inv(material.rho)+s*ε0*material.eps_r)-gamma0
            @test recovered ≈ equivalent^2 rtol=1e-12
            functor=E.EarthAdmittance.Formula(:Xue2021)(
                [Inf,material.rho],[ε0,ε0*material.eps_r],[μ0,μ0],s,nothing)
            hi,hj=pair.heights; H=hi+hj; x=pair.separation
            D=hypot(H,x); d=hypot(hi-hj,x)
            correction=quadgk(t->exp(-H*t)*cos(x*t)/
                (equivalent^2/gamma0*t+sqrt(t^2+equivalent^2)),
                0.0,Inf;rtol=1e-10)[1]
            expected=(log(D/d)+2correction)/(2π*ε0)
            @test functor(kind,pair) ≈ expected rtol=1e-9
            @test isfinite(functor(kind,pair))
            reversepair=E.EarthPair(pair.column,pair.row,reverse(pair.heights),x,(1,1))
            @test functor(kind,reversepair) ≈ expected rtol=1e-9
        end
    end
    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    design=build(CableDesign,"EHEM-potential-wire",
        Group(:core,Region(:core,Disk(0.01),copper)))
    system=build(LineCableSystem,[design,design],[(0.0,10.0),(3.0,15.0)];
        connections=[Dict("core"=>1),Dict("core"=>2)])
    problem=LineParametersProblem(system;earth_props=model,frequencies=[50.0,1e5])
    formulation=Formulation(earth_impedance=:Wise1934,earth_admittance=:Xue2021,
        equivalent_earth=formula(:Xue2021;order=:after),
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    result=compute(problem,formulation;options=(trace=true,))
    @test size(result.Y)==(2,2,2)
    @test all(isfinite,result.Y)
    @test details(result).trace.Pg[:,:,1] ≈ transpose(details(result).trace.Pg[:,:,1])
end
