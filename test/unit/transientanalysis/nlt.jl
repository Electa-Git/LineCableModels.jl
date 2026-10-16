@testitem "TransientAnalysis / NLT exact lossless line" tags=[:unit, :engine] setup=[TestFixtures] begin
    using Statistics: mean
    import LineCableModels.TransientAnalysis as TA
    import LineCableModels.Engine: SeriesImpedance, ShuntAdmittance, PhaseDomain

    design = build(CableDesign, "overhead", Stack(
        Group(:core, Region(:metal, Disk(0.01), TestFixtures.conductor_material()))))
    system = build(LineCableSystem, [design], [Pose2(0.0, 10.0)];
        connections=[Dict(:core=>1)], line_length=1000.0)
    line = LineParametersProblem(system; earth_props=homogeneous(rho=100.0))
    window = (tmax=4e-5, step=1e-7)
    terminal = "cable:1:core"
    selected = TransientFormulation(:nlt; options=(N=2048,)).formula
    N = selected.options.data.N
    T = selected.options.data.Tfactor * window.tmax
    c = log(N^2) / T
    f = [(2k - 1) / (2T) for k in 1:N]
    s = complex.(c, 2π .* f)
    L, C = 1e-6, 1e-10
    scan = LineParameters(PhaseDomain(),
        SeriesImpedance{ComplexF64,:pul}(reshape(s .* L, 1, 1, N)),
        ShuntAdmittance{ComplexF64,:pul}(reshape(s .* C, 1, 1, N)), f,
        ComputationDetails(; coordinates=[terminal]))
    arrival = 1000sqrt(L * C)
    for (receiving, expected) in (((), 2.0), (Matched(), 1.0))
        problem = TransientAnalysisProblem(line; sending=(Step(terminal, 1.0),), receiving, window)
        workspace = TA.TransientWorkspace(scan, problem, selected)
        t, voltages, currents = TA.respond!(selected, workspace, selected.parameters, selected.options)
        @test workspace.input.ℓ == 1000
        @test workspace.plan.coordinates == [terminal]
        @test size(voltages) == size(currents) == (2, 1, length(t))
        @test maximum(abs, voltages[2, 1, t .< arrival - 2window.step]) < 1e-3
        plateau = (1.2arrival .<= t) .& (t .<= 2.8arrival)
        @test mean(voltages[2, 1, plateau]) ≈ expected rtol=1e-2
    end
    problem = TransientAnalysisProblem(line; sending=(Step(terminal, 1.0),), receiving=(), window)
    reduced = LineParameters(PhaseDomain(), scan.Z, scan.Y, f,
        ComputationDetails(; coordinates=["removed"]))
    @test_throws ArgumentError("terminal $terminal: the line formulation removed it") TA.TransientWorkspace(
        reduced, problem, selected)
    @test_throws ArgumentError TransientAnalysisProblem(line;
        sending=(Step("unknown", 1.0),), receiving=(), window)
    @test_throws DomainError TransientFormulation(:nlt; options=(N=1,))
    for Tfactor in (1.0, Inf, NaN)
        @test_throws DomainError TransientFormulation(:nlt; options=(Tfactor=Tfactor,))
    end
end
