@testitem "TransientAnalysis / raw observations and text display" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels.TransientAnalysis as TA

    t = [0.0, 0.1, 0.2]
    voltages = reshape(collect(1.0:12.0), 2, 2, 3)
    currents = -voltages
    response = TransientResponse(t, voltages, currents,
        ["cable:1:core", "cable:2:core"], ComputationDetails())
    for (selector, values) in ((voltage, voltages), (current, currents)),
            (index, label) in ((1, "sending"), (2, "receiving"))
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
