@testitem "ObservedResult / operand eligibility precedes RMS" tags=[:unit] begin
    using LineCableModels.Engine: compare
    using LineCableModels.Grammar: observation_resolution
    using LinearAlgebra: norm
    for T in (Float32,Float64,BigFloat)
        cutoff=T(1//10^10)
        frequency=T[1,10,100]
        resistance=reshape(T[-2cutoff,-cutoff,zero(T)],1,1,3)
        classified=observation_resolution(resistance,R;atol=(R=cutoff,),frequencies=frequency)
        @test vec(classified.unresolved)==[false,true,true]
        @test all(classified.available)
        @test all(observation_resolution(zero.(resistance),R;atol=0).unresolved)
        a=LineParameters(PhaseDomain,complex.(resistance,one(T)),fill(complex(one(T),one(T)),1,1,3),frequency)
        b=LineParameters(PhaseDomain,complex.(2resistance,one(T)),a.Y,frequency)
        for normalization in (:reference_rms,:pointwise)
            result=compare(a,b,R;atol=(R=cutoff,),normalization)
            @test ismissing(only(result.absolute))
            @test ismissing(only(result.relative))
        end
        # A valid error can be smaller than the physical operand floor.
        original=fill(T(1e-3),1,1,3)
        changed=original .+ T(1e-11)
        if any(changed .!= original)
            result=compare(complex.(original,zero(T)),complex.(changed,zero(T)),R;
                frequencies=frequency,result_basis=:pul,atol=(R=cutoff,))
            @test 0 < only(result.absolute) < cutoff
            @test only(result.absolute) ≈ norm(vec(changed-original))/sqrt(T(3))
        end
        @test only(compare(a,a,R;band=(1,1),atol=(R=cutoff,)).absolute)==0
        @test only(compare(a,a,R;band=(1,1),atol=(R=cutoff,)).relative)==0
    end
    # Components, not a magnitude approximation, decide complex zero.
    z=reshape(ComplexF64[8e-11+8e-11im, 1.1e-10+0im],1,1,2)
    resolution=observation_resolution(z,Z;frequencies=[1.,2.],atol=(R=1e-10,X=1e-10))
    @test vec(resolution.unresolved)==[true,false]
    @test_throws ArgumentError observation_resolution(z,Z;atol=(Z=1e-10,))
    @test_throws ArgumentError observation_resolution(z,Z;atol=(X=1e-10,L=1e-15,))
    a=LineParameters(PhaseDomain,fill(1.0+1im,1,1,2),fill(1.0+1im,1,1,2),[0.,1.])
    @test ismissing(only(compare(a,a,L).absolute))
    @test only(compare(a,a,R).absolute)==0
    a.Z.values[1,1,1]=complex(NaN,0)
    @test ismissing(only(compare(a,a,R).absolute))
    @test only(compare(a,a,R;band=(1.,1.)).absolute)==0
end

@testitem "ObservedResult / declared inputs survive completion without tracing" tags=[:unit] setup=[TestFixtures] begin
    using LineCableModels.Grammar: observation_gridpoint
    problem=TestFixtures.line_parameters_problem()
    scalar=compute(problem)
    retained=observation_gridpoint(scalar)
    @test retained.inputs.temperature==20.0
    @test retained.inputs.earth_props.layers[2].rho==100.0
    @test retained.inputs.system.line_length==600.0
    @test retained.inputs.system.designs[1].origin.items[1].item.items[1].primitive.r==0.005
    @test retained.id.problem_index==retained.id.formulation_index==1
    @test retained.formulations==scalar.details.data.formulations
    calls=Ref(0)
    space=Gridspace{LineParametersProblem}(rho -> begin
        calls[]+=1
        LineParametersProblem(problem.system;earth_props=EarthModel(rho),frequencies=[50.])
    end,(Grid([100.,200.]),))
    results=compute(ParametricProblem(space),Combinatorial(Formulation()))
    @test calls[]==2
    descriptions=observation_gridpoint.(collect(results))
    @test calls[]==2
    @test [d.id.problem_index for d in descriptions]==[1,2]
    @test [d.inputs.earth_props.layers[2].rho for d in descriptions]==[100.,200.]
    @test all(d -> d.inputs.system.designs[1].origin.items[1].item.items[1].primitive.r==0.005,descriptions)
    descriptions[1].inputs.system.connection_order[1]=-100
    @test results[1].details.data.inputs.system.connection_order[1] != -100
    @test observation_gridpoint(results[2]).inputs.earth_props.layers[2].rho==200.
    constants=compute(CableConstantsProblem(problem.system.designs[1]))
    @test observation_gridpoint(constants).inputs.design.origin.items[1].item.items[1].primitive.r==0.005
end

@testitem "ObservedResult / nominal recentering retains uncertainty dependencies" tags=[:unit] begin
    using Measurements, Calculus
    using LineCableModels.Grammar: observation_resolution
    shared=Measurements.measurement(5e-13,1e-20)
    x=reshape([shared,2shared],1,1,2)
    resolution=observation_resolution(x,G;frequencies=[1.,2.])
    shifted=LineCableModels.Grammar._resolved_observation(x,resolution.unresolved,resolution.available,Val(false))
    @test Measurements.value.(shifted)==zeros(1,1,2)
    @test Measurements.uncertainty.(shifted)==Measurements.uncertainty.(x)
    @test Measurements.uncertainty(shifted[2]-2shifted[1])==0
    @test Measurements.uncertainty(shifted[1]-shared)==0
    # Native Measurements.isfinite checks the nominal value, not its spread.
    undefined=abs(complex(Measurements.measurement(0.,1.),Measurements.measurement(0.,2.)))
    @test !LineCableModels.Engine._resolution_available(undefined)
end
