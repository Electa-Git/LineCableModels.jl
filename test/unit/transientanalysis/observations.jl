@testitem "TransientAnalysis / raw observations and text display" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels.TransientAnalysis as TA

    t = [0.0, 0.1, 0.2]
    voltages = reshape(collect(1.0:12.0), 2, 2, 3)
    currents = -voltages
    response = TransientResponse(t, voltages, currents,
        ["cable:1:core", "cable:2:core"], ComputationDetails())
    for (selector, values) in ((voltage, voltages), (current, currents)),
            (index, label) in ((1, "sending"), (2, "receiving"))
        @test observe(response, selector) === values
        @test observe(response, selector, index, 1, :) == values[index, 1, :]
        @test observe(response, selector, label, "cable:2:core", :) == values[index, 2, :]
        @test observe(response, selector, index, "cable:1:core", :) == values[index, 1, :]
        @test observe(response, selector, label, 2, :) == values[index, 2, :]
    end
    @test_throws ArgumentError observe(response, voltage, "unknown", 1, :)
    @test_throws ArgumentError observe(response, current, 1, "unknown", :)
    @test_throws BoundsError observe(response, voltage, 3, 1, :)
    @test_throws BoundsError observe(response, current, 1, 3, :)

    problem = TransientAnalysisProblem(TestFixtures.three_bare_wires_problem();
        sending=(Step("cable:1:core", 1),), receiving=Matched(),
        window=(tmax=0.2, step=0.1))
    formulation = TransientFormulation()
    workspace = TA.TransientWorkspace((;), (coordinates=response.coordinates, window=problem.window), (;))
    for (value, name) in ((problem, "TransientAnalysisProblem"),
            (formulation, "TransientFormulation"), (workspace, "TransientWorkspace"),
            (response, "TransientResponse"))
        @test occursin(name, sprint(show, value))
        @test !isempty(sprint(summary, value))
        @test !isempty(sprint(show, MIME"text/plain"(), value))
        @test sprint(show, MIME"text/plain"(), value; context=:compact=>true) == sprint(show, value)
    end
    @test occursin(":nlt", sprint(show, formulation))
    @test occursin("voltage", sprint(show, MIME"text/plain"(), response))
    @test occursin("current", sprint(show, MIME"text/plain"(), response))
end

@testitem "TransientAnalysis / retained observations, reports and archives" tags=[:unit, :engine] begin
    using DataFrames
    import LineCableModels.Commons as G
    import LineCableModels.Units as U

    t=[0.0,1e-6,2e-6]
    voltages=reshape(collect(1.0:12.0),2,2,3)
    currents=-voltages
    response=TransientResponse(t,voltages,currents,
        ["cable:1:core","cable:2:core"],ComputationDetails())
    observed=ObservedResult(response)
    @test G.observables(typeof(response))==(voltage,current)
    @test G.basis(response)===:total
    @test length(observed.quantities)==2
    @test observe(observed,voltage)==voltages
    @test observe(observed,current)==currents
    @test observed.quantities[1].values!==voltages
    @test observed.quantities[1].coordinates.times==t
    @test observed.quantities[1].coordinates.times!==t
    @test observed.quantities[1].coordinates.column_labels!==response.coordinates
    @test all(q -> q.thresholds===q.available===q.engineering_zero===q.missing_reason===nothing &&
        !q.clipped && q.family===:transient && q.statistic===:value,observed.quantities)
    @test ObservedResult(response;clip=false).quantities==observed.quantities
    @test_throws ArgumentError ObservedResult(response;atol=0)
    @test_throws ArgumentError ObservedResult(response;frequencies=[1.,2.,3.])
    @test_throws DimensionMismatch ObservedResult(response,((voltage,1,1),))
    @test_throws BoundsError ObservedResult(response,((current,3,:,:),))

    for source in (response,observed)
        micro=ObservedResult(source;time_unit=:micro)
        @test micro.quantities[1].coordinates.times==t.*1e6
        @test micro.quantities[1].coordinates.time_unit==U.units(:micro,:second)
        @test observe(micro,voltage)==voltages
        @test G.observation_product(micro,voltage).coordinates.times==t.*1e6
        @test G.observation_product(micro,voltage).coordinates.time_unit==U.units(:micro,:second)
        @test ObservedResult(micro).quantities==micro.quantities
        @test G.observation_product(micro,voltage;time_unit=:micro).coordinates.times==t.*1e6
        @test ObservedResult(micro;time_unit=:base).quantities==observed.quantities
    end
    for request in ((voltage,2,1,2),(voltage,1:2,2,2:3),(current,:,1:2,:))
        indexed=ObservedResult(response,(request,))
        retained=G.observation_product(observed,request)
        @test observe(indexed,request)==retained.values
        @test retained.values==observe(response,first(request))[request[2:end]...]
        @test indexed.quantities[1].coordinates==retained.coordinates
    end
    selection=(voltage,[2,1],[2,1],[3,1])
    subset=ObservedResult(response,(selection,))
    original=(voltage,1,2,1)
    @test observe(subset,original)==voltages[1,2,1]
    @test G.observation_product(subset,original).coordinates.samples==[1]
    @test G.observation_product(subset,original).coordinates.times==[t[1]]
    @test_throws ArgumentError observe(subset,voltage,1,1,2)
    displayed=ObservedResult(response,(voltage,);quantity_units=:kilo)
    @test observe(displayed,voltage)==voltages.*1e-3

    artifact=report(TableReportDefinition((voltage,current)),observed)
    table=artifact[voltage]
    @test names(table)==["time","sending cable:1:core","sending cable:2:core",
        "receiving cable:1:core","receiving cable:2:core"]
    @test table.time==t
    @test table[!,2]==voltages[1,1,:]
    @test metadata(table,"observation_columns").time==
        (quantity=U.Quantity{:time}(),unit=U.units(:base,:second))
    @test occursin("sending cable:1:core",sprint(show,MIME"text/plain"(),artifact;context=:limit=>false))
    @test occursin("time",sprint(show,MIME"text/html"(),artifact))
    for reference in (observed,response)
        converted=report(response,(voltage,);reference,time_unit=:micro)
        @test converted.observed.quantities[1].coordinates.times==t.*1e6
        @test converted.reference.quantities[1].coordinates.times==t.*1e6
    end
    mktempdir() do directory
        for extension in ("json","jls")
            path=joinpath(directory,"transient.$extension")
            save(observed,path)
            restored=import_data(:observed,path)
            @test restored.gridpoint==observed.gridpoint
            @test restored.quantities==observed.quantities
            @test restored.errors==observed.errors
            @test restored.timings==observed.timings
            @test observe(restored,voltage)==voltages
            @test observe(restored,current)==currents
            @test report(restored)[voltage]==table
        end
    end
end
