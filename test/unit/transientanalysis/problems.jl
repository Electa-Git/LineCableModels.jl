@testitem "TransientAnalysis / problem and formulation grammar" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels: validate
    import LineCableModels.TransientAnalysis as TA

    line = TestFixtures.three_bare_wires_problem()
    window = (tmax=1e-3, step=1e-5)
    sending = (Step("cable:1:core", 1),)
    problem = TransientAnalysisProblem(line; sending, receiving=(), window)
    @test problem.line === line
    @test problem.sending === sending
    @test problem.receiving === ()
    @test problem.window === window
    @test validate(problem) === problem
    @test_throws ArgumentError TransientAnalysisProblem(line;
        sending=(Step("unknown", 1),), receiving=(), window)
    @test_throws ArgumentError TransientAnalysisProblem(line; sending,
        receiving=(TerminalImpedance("unknown"),), window)
    @test_throws ArgumentError TransientAnalysisProblem(line;
        sending=(sending..., TerminalImpedance("cable:1:core")), receiving=(), window)
    @test_throws ArgumentError TransientAnalysisProblem(line; sending,
        receiving=(Step("cable:2:core", 1), TerminalImpedance("cable:2:core")), window)
    passive = TransientAnalysisProblem(line;
        sending=(TerminalImpedance("cable:1:core"),), receiving=(), window)
    @test validate(passive) === passive
    matched = TransientAnalysisProblem(line;
        sending=Matched(), receiving=Matched(), window)
    @test validate(matched) === matched
    open = TransientAnalysisProblem(line; sending=(), receiving=(), window)
    @test validate(open) === open
    for bad in ((tmax=0.0, step=1e-5), (tmax=-1.0, step=1e-5),
            (tmax=Inf, step=1e-5), (tmax=NaN, step=1e-5),
            (tmax=1e-3, step=0.0), (tmax=1e-3, step=-1e-5),
            (tmax=1e-3, step=1e-3), (tmax=1e-3, step=1.0),
            (tmax=1e-3, step=Inf), (tmax=1e-3, step=NaN))
        @test_throws DomainError TransientAnalysisProblem(line; sending, receiving=(), window=bad)
    end
    for (left, right) in ((sending, Matched()), (Matched(), sending))
        matched = TransientAnalysisProblem(line; sending=left, receiving=right, window)
        @test validate(matched) === matched
    end
    space = TransientAnalysisProblem(line;
        sending=(sending..., TerminalImpedance("cable:2:core"; R=Grid((1, 2)))),
        receiving=Matched(), window)
    @test space isa Gridspace{TransientAnalysisProblem}
    @test [point.sending[2].R for point in space] == [1, 2]
    @test TransientAnalysisProblem(Grid((line, line)); sending, receiving=(), window) isa Gridspace
    @test TransientAnalysisProblem(line; sending, receiving=(), window=Grid((window, window))) isa Gridspace

    default = TransientFormulation()
    @test formula_id(default.formula) === :nlt
    controls = formulation_options(default.formula).data
    @test controls.line isa LineParametersFormulation
    @test (controls.N, controls.Tfactor, controls.taper) == (4096, 1.25, :hanning)
    @test TA.formulas(TA.Formula) == (:nlt, :td, :default)
    for taper in (:hanning, :blackman, :lanczos)
        selected = TransientFormulation(:nlt; options=(N=128, Tfactor=2.0, taper=taper))
        @test formulation_options(selected.formula).data.taper === taper
    end
    @test isempty(formulation_options(TransientFormulation(:td).formula).data)
    @test TransientFormulation(Grid((:nlt, :td))) isa Gridspace{TransientFormulation}
    @test [formula_id(item.formula) for item in TransientFormulation(Grid((:nlt, :td)))] == [:nlt, :td]
    @test_throws ArgumentError TransientFormulation(:nlt; options=(unknown=1,))
    @test_throws ArgumentError TransientFormulation(:td; options=(N=128,))
    @test_throws ArgumentError TransientFormulation(:UnregisteredIdentity)
    @test_throws TypeError TransientFormulation(:nlt; options=(N=1.0,))
    @test_throws TypeError TransientFormulation(:nlt; options=(line=nothing,))
    @test_throws ArgumentError TransientFormulation(:nlt; options=(taper=:unknown,))
    execution = computation_options(TransientFormulation, ComputationOptions())
    @test execution.data == (verbosity=(default=0,), timing=false)
    @test_throws ArgumentError computation_options(TransientFormulation, ComputationOptions(; unknown=true))
    @test_throws TypeError computation_options(TransientFormulation, ComputationOptions(; timing=1))
    @test_throws ArgumentError computation_options(TransientFormulation, ComputationOptions(; verbosity=(default=3,)))
    @test TA.initialize_buffers(default.formula, Float64, nothing, nothing, (;)) === (;)
    @test validate(problem, TransientFormulation(:nlt).formula) === problem
    for identifier in (:td,)
        selected = TransientFormulation(identifier)
        @test_throws ArgumentError("not implemented yet") TA.respond!(selected.formula, nothing, (;), selected.formula.options)
        message = try
            validate(problem, selected.formula)
            ""
        catch error
            sprint(showerror, error)
        end
        @test occursin(string(identifier), message)
        @test occursin(string(typeof(line)), message)
        @test_throws ArgumentError compute(problem, selected)
    end
end
