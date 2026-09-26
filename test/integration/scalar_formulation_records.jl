@testitem "ReportBuilder / scalar results retain physical and native formula choices" tags=[:integration] setup=[FormulaFixtures] begin
    using LinearAlgebra, Serialization
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    copper=Material(MaterialsLibrary(add_defaults = true), :copper)
    design=build(CableDesign, "scalar-records", terminal(:core, solid(copper, Disk(0.0425))))
    system=build(LineCableSystem, [design, design], [Pose2(0, -1), Pose2(1, -1)];
        connections = [Dict(:core=>1), Dict(:core=>2)])
    problem=LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = [1e5])
    # Material coefficients are selectable; receiving-layer voltage references
    # are determined by the physical geometry, not formulation parameters.
    coefficients=(0.1, 1.0)
    choices=[Formulation(earth_properties = formula(:portela1999; parameters = (; beta)))
             for beta in coefficients]
    results=compute.(Ref(problem), choices)
    @test norm(Y(results[1])-Y(results[2]))/norm(Y(results[1])) > 0.1
    @test typeof(results[1]) === typeof(results[2])
    for (result, selection, beta) in zip(results, choices, coefficients)
        record=details(result).data.formulations
        @test record.requested == NamedTuple(selection).requested
        @test record.methods == NamedTuple(selection).methods
        @test record.requested.earth_properties.parameters.beta === beta
        @test record.methods.earth_properties.parameters.beta === beta
    end
    artifact=report(BenchmarkTableDefinition(), (
        reference = results[1], candidate = results[2]))
    @test artifact.reference.gridpoint.formulations !=
          artifact.observed.gridpoint.formulations
    @test artifact.tables.formulations.label[1] != artifact.tables.formulations.label[2]
    @test all(!isempty, artifact.tables.formulations.label)
    io=IOBuffer();
    serialize(io, results);
    seekstart(io)
    restored=deserialize(io)
    # Unavailable cells must stay unavailable after serialization; == propagates
    # missing instead of comparing the retained availability mask.
    @test isequal(
        report(BenchmarkTableDefinition(), (
            reference = restored[1], candidate = restored[2])).tables.formulations,
        artifact.tables.formulations)
    custom=FormulaFixtures.DispersiveEarth()
    selected=Formulation(earth_properties = custom)
    changed=compute(problem, selected)
    @test !isempty(custom.seen)
    @test details(changed).data.formulations.requested.earth_properties==NamedTuple(custom)
    @test details(changed).data.formulations.methods.earth_properties==NamedTuple(custom)
    @test details(changed).data.formulations.methods.earth_properties.identifier === :DispersiveEarth
    @test typeof(changed) === typeof(results[1])
    retained=LineCableModels.ImportExport.deserialize_value(
        Val(:formulation), details(changed).data.formulations)
    @test formula_id(retained, nothing)==formula_id(selected, nothing)
    grid=Formulation(
        earth_properties = Grid([c.definitions.earth_properties for c in choices]);
        combine = :zip)
    batch=compute(problem, grid)
    @test isconcretetype(eltype(batch))
    @test length(batch)==2
    for i in 1:2
        @test Z(batch[i]) == Z(results[i])
        @test Y(batch[i]) == Y(results[i])
        @test details(batch[i]).data.formulations == details(results[i]).data.formulations
    end
end
