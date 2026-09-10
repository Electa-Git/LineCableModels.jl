@testitem "ReportBuilder / scalar results retain physical formula choices and hooks" tags=[:integration] begin
    using LinearAlgebra, Serialization
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    copper=Material(MaterialsLibrary(add_defaults=true),:copper)
    design=build(CableDesign,"scalar-records",terminal(:core,solid(copper,Disk(0.0425))))
    system=build(LineCableSystem,[design,design],[Pose2(0,-1),Pose2(1,-1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem=LineParametersProblem(system;earth_props=homogeneous(rho=0.1),frequencies=[1e3])
    choices=[Formulation(earth_impedance=formula(:default;parameters=(;reference)),
        earth_admittance=formula(:default;parameters=(;reference))) for reference in (:deep,:interface)]
    results=compute.(Ref(problem),choices)
    @test norm(Y(results[1])-Y(results[2]))/norm(Y(results[1])) > 0.1
    @test typeof(results[1]) === typeof(results[2])
    for (result,selection,reference) in zip(results,choices,(:deep,:interface))
        record=details(result).formulations
        @test record.requested == NamedTuple(selection).requested
        @test record.methods == NamedTuple(selection).methods
        @test record.requested.earth_admittance.parameters.reference === reference
        @test record.methods.earth_admittance.parameters.reference === reference
    end
    artifact=report(BenchmarkTableDefinition(),(reference=results[1],candidate=results[2]))
    @test artifact.table.formulations.record[1] != artifact.table.formulations.record[2]
    @test occursin("deep",artifact.table.formulations.label[1])
    @test occursin("interface",artifact.table.formulations.label[2])
    io=IOBuffer();serialize(io,results);seekstart(io)
    restored=deserialize(io)
    @test report(BenchmarkTableDefinition(),(reference=restored[1],candidate=restored[2])).table.formulations == artifact.table.formulations
    calls=Ref(0)
    hook=(s,m,l)->(calls[]+=1;zero(s))
    hooked=compute(problem,Formulation(earth_impedance=formula(:default;hooks=(Γ=hook,))))
    @test calls[] > 0
    @test details(hooked).formulations.requested.earth_impedance.hooks.Γ === hook
    @test details(hooked).formulations.methods.earth_impedance.hooks.Γ === hook
    @test typeof(hooked) === typeof(results[1])
    grid=Formulation(earth_impedance=Grid([c.definitions.earth_impedance for c in choices]),
        earth_admittance=Grid([c.definitions.earth_admittance for c in choices]);combine=:zip)
    batch=compute(problem,grid)
    @test isconcretetype(eltype(batch))
    @test length(batch)==2
    for i in 1:2
        @test Z(batch[i]) == Z(results[i])
        @test Y(batch[i]) == Y(results[i])
        @test details(batch[i]).formulations == details(results[i]).formulations
    end
end
