@testitem "TransientAnalysis / NLT three bare wires" tags=[:integration, :engine] setup=[TestFixtures] begin
    line = TestFixtures.three_bare_wires_problem(; line_length=1000.0)
    sending = (Step("cable:1:core", 1.0),
        TerminalImpedance("cable:2:core"; R=1.0),
        TerminalImpedance("cable:3:core"; R=1.0))
    receiving = (TerminalImpedance("cable:2:core"; R=1.0),
        TerminalImpedance("cable:3:core"; R=1.0))
    problem = TransientAnalysisProblem(line; sending, receiving,
        window=(tmax=2e-4, step=1e-6))
    formulation = TransientFormulation(:nlt; options=(N=256,
        line=Formulation(earth_impedance=:unified, earth_admittance=:unified)))
    response = compute(problem, formulation; options=(timing=true,))
    @test all(isfinite, response.voltage)
    @test all(isfinite, response.current)
    @test response.coordinates == ["cable:$i:core" for i in 1:3]
    @test maximum(abs, response.voltage[1, 1, response.t .> 5problem.window.step] .- 1) < 1e-2
    largest = maximum(abs, response.current)
    @test all(abs.(response.current[:, 2:3, :] + response.voltage[:, 2:3, :]) .<= 1e-9largest)
    retained = details(response).data
    @test retained.nlt.N == 256
    @test retained.nlt.T == 1.25problem.window.tmax
    @test retained.nlt.c == log(256^2) / retained.nlt.T
    @test retained.nlt.Tfactor == 1.25
    @test retained.nlt.taper === :hanning
    @test retained.nlt.aliasing ≈ exp(-2retained.nlt.c * retained.nlt.T) /
        (1 - exp(-2retained.nlt.c * retained.nlt.T))
    @test retained.line isa ComputationDetails
    @test retained.timing.wall_seconds >= 0
end
