@testitem "ReportBuilder / scalar results retain physical and native formula choices" tags=[:integration, :importexport, :slow] setup=[FormulaFixtures] begin
    using LinearAlgebra, Serialization, JSON3
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    copper=Material(MaterialsLibrary(add_defaults = true), :copper)
    design=build(CableDesign, "scalar-records", terminal(:core, solid(copper, Disk(0.0425))))
    system=build(LineCableSystem, [design, design], [Pose2(0, -1), Pose2(1, -1)];
        connections = [Dict(:core=>1), Dict(:core=>2)])
    problem=LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = [1e5])
    # Material coefficients are selectable. Receiving-layer voltage references
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
        reference = results[1], result = results[2]))
    @test artifact.reference.gridpoint.formulations !=
          artifact.observed.gridpoint.formulations
    @test artifact.tables.formulations.label[1] != artifact.tables.formulations.label[2]
    @test all(!isempty, artifact.tables.formulations.label)
    io=IOBuffer();
    serialize(io, results);
    seekstart(io)
    restored=deserialize(io)
    # Unavailable cells must remain unavailable after serialization. == propagates
    # missing instead of comparing the retained availability mask.
    @test isequal(
        report(BenchmarkTableDefinition(), (
            reference = restored[1], result = restored[2])).tables.formulations,
        artifact.tables.formulations)
    custom=FormulaFixtures.DispersiveEarth()
    selected=Formulation(earth_properties = custom)
    changed=compute(problem, selected)
    @test !isempty(custom.seen)
    @test details(changed).data.formulations.requested.earth_properties==NamedTuple(custom)
    @test details(changed).data.formulations.methods.earth_properties==NamedTuple(custom)
    @test details(changed).data.formulations.methods.earth_properties.identifier === :DispersiveEarth
    @test typeof(changed) === typeof(results[1])
    # A saved and reloaded result has the details it was saved with, labels included, for a
    # built-in and for a custom formula. Loading builds no formulation.
    IE=LineCableModels.ImportExport
    for saved in (results[1], changed)
        loaded=IE.deserialize_value(JSON3.read(JSON3.write(IE.serialize_value(saved)), Dict{String,Any}))
        @test isequal(details(loaded).data, details(saved).data)
        @test typeof(details(loaded).data) === typeof(details(saved).data)
        @test details(loaded).data.formulation_fields == details(saved).data.formulation_fields
        @test Z(loaded) == Z(saved) && Y(loaded) == Y(saved)
    end
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
