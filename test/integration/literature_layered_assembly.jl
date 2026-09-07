@testitem "Engine / literature assimilation / layered source terms in full matrices" begin
    using LineCableModels, LinearAlgebra
    E=LineCableModels.Engine; EP=LineCableModels.EarthProps
    copper=Material(:conductor,1.7241e-8,1.,1.,20.,0.)
    wire=build(CableDesign,"layered wire",Group(:core,Region(:core,Disk(.01),copper)))
    system=build(LineCableSystem,[wire,wire],[(0.,10.),(3.,15.)];
        connections=[Dict("core"=>1),Dict("core"=>2)])
    μ0,ε0=4π*1e-7,8.8541878128e-12
    for (id,n) in ((:Wedepohl1966,2),(:Moghram1998,3),(:Lee2014,3))
        rho=[100.,500.,50.][1:n]; er=[10.,6.,20.][1:n]
        thickness=[2.,3.][1:n-1]
        earth=build(EarthModel,Tuple(EP.EarthLayer(rho[k],er[k],1.,
            k==n ? Inf : thickness[k]) for k in 1:n))
        problem=LineParametersProblem(system;earth_props=earth,frequencies=[50.,1e5])
        methods=Formulation(earth_impedance=id,earth_admittance=:Ametani2021,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        result=compute(problem,methods;options=(trace=true,))
        trace=details(result).trace
        for (k,freq) in enumerate(problem.frequencies)
            s=complex(0.,2π*freq)
            leaf=E.EarthImpedance.Formula(id)([Inf;rho],ε0.*[1.;er],
                fill(μ0,n+1),s,nothing,nothing,[Inf;thickness;Inf])
            expected=zeros(ComplexF64,2,2)
            for i in 1:2,j in 1:2
                pair=E.EarthPair(i,j,([10.,15.][i],[10.,15.][j]),
                    i==j ? .01 : 3.,(1,1))
                expected[i,j]=leaf(i==j ? Val(:self) : Val(:mutual),pair)
            end
            @test trace.Zg[:,:,k] ≈ expected rtol=2e-8
            @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k])
            @test result.Y[:,:,k] ≈ s*inv(trace.P[:,:,k]) rtol=2e-10
            @test all(isfinite,result.Z[:,:,k])
            @test minimum(eigvals(Symmetric(real.(result.Z[:,:,k]))))>0
        end
    end
end
