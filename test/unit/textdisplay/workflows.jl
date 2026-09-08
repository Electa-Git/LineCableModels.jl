@testitem "TextDisplay / calculation intent and numerical results" tags=[:unit] setup=[
    TestFixtures,
] begin
    const EN=LineCableModels.Engine
    design=TestFixtures.mv_cable_design()
    problem=CableConstantsProblem(design; temperature = 60.0, frequency = 60.0)
    line_problem=LineParametersProblem(TestFixtures.three_phase_system();
        earth_props = homogeneous(rho = 100.0), frequencies = [0.1, 50.0, 1e6])
    formulation=Formulation(
        earth_impedance = :default,
        insulation_admittance = :Ametani2004,
        semicon_admittance = :default
    )
    cable_formulation=CableConstantsFormulation(insulation_admittance = :Ametani2004)
    parameters=TestFixtures.two_conductor_results()
    modal=compute(ModalTransformationProblem(parameters), ModalTransformationFormulation())
    benchmark=EN.compare(parameters, parameters)
    constants=CableConstants(1e-4, 2e-7, 3e-10, 4e-12)
    backend=LineCableModelsFEM(fem_options = (mesh_policy = :remesh,))
    blueprints=EN.CableBlueprint{eltype(line_problem)}[EN.flatten(LineCableModelsCoaxial(),
                                                           design, eltype(line_problem))
                                                       for design in line_problem.system.designs]
    workspaces=map((false, true)) do trace
        execution=LineCableModels.Grammar.computation_options(LineCableModelsCoaxial, (;
            trace))
        EN.LineParametersWorkspace(line_problem, formulation, execution, blueprints)
    end
    for workspace in workspaces
        # Scratch storage is intentionally uninitialized until compute. Give
        # the two primitive buffers sentinels before checking display purity.
        fill!(workspace.buffers.Zbuffer, 123+im)
        fill!(workspace.buffers.Pbuffer, -456+2im)
    end
    objects=(
        problem, line_problem, formulation, cable_formulation, constants,
        LineCableModelsCoaxial(), backend, backend.execution,
        EN.PhaseDomain(), modal.domain, parameters.Z, parameters.Y,
        parameters, modal, benchmark, benchmark.Z, workspaces...
    )

    # REPL inspection must identify the requested computation, not dump matrix
    # contents, solve a problem, or require the optional backend to be loaded.
    before=(copy(parent(parameters.Z)), copy(parent(parameters.Y)))
    for object in objects
        compact=sprint(show, object)
        @test !occursin('\n', compact)
        @test !isempty(sprint(summary, object))
        @test sprint(show, MIME"text/plain"(), object; context = :compact=>true) == compact
        for display_size in ((40, 120), (5, 48))
            context=IOContext(IOBuffer(), :limit=>true, :displaysize=>display_size)
            shown=sprint(show, MIME"text/plain"(), object; context)
            @test !endswith(shown, '\n')
            @test length(split(shown, '\n')) <= first(display_size)
            @test all(textwidth(line) <= last(display_size) for line in split(shown, '\n'))
        end
    end
    @test parent(parameters.Z) == before[1]
    @test parent(parameters.Y) == before[2]
    for (index, workspace) in enumerate(workspaces)
        @test all(==(123 + im), workspace.buffers.Zbuffer)
        @test all(==(-456 + 2im), workspace.buffers.Pbuffer)
        @test occursin("frequencies=3", sprint(show, workspace))
        @test occursin(index == 1 ? "disabled" : "enabled",
            sprint(show, MIME"text/plain"(), workspace))
    end
    shown=sprint(show, MIME"text/plain"(), problem)
    @test occursin(design.cable_id, shown)
    @test occursin("60 °C", shown)
    @test occursin("60 Hz", shown)
    shown=sprint(show, MIME"text/plain"(), line_problem)
    @test occursin("3 points", shown)
    @test occursin("100 mHz", shown)
    @test occursin("1 MHz", shown)
    shown=sprint(show, MIME"text/plain"(), formulation;
        context = IOContext(IOBuffer(), :displaysize=>(40, 160)))
    @test occursin("default", shown)
    @test occursin("Ametani2004", shown)
    @test occursin("semicon_admittance", shown)
    @test occursin("default", shown)
    @test occursin("Ametani2004", sprint(show, MIME"text/plain"(), cable_formulation))
    @test occursin("Ω/m", sprint(show, MIME"text/plain"(), parameters.Z))
    @test occursin("S/m", sprint(show, MIME"text/plain"(), parameters.Y))
    @test occursin("RMS errors", sprint(show, MIME"text/plain"(), benchmark))
    @test occursin("relative", sprint(show, MIME"text/plain"(), benchmark.Z))
    @test occursin("G=", sprint(show, MIME"text/plain"(), constants))
    @test occursin("remesh", sprint(show, backend.execution))
    @test occursin("modal domain", sprint(show, modal))
    @test occursin("Modal domain", sprint(summary, modal.domain))
end

@testitem "TextDisplay / wire estimates expose feasibility and physical dimensions" tags=[:unit] begin
    estimates = (
        make_stranded(1000.0),
        make_stranded(1.0e12; nmin = 40, nmax = 40),
        make_screened(35.0, 60.0; coverage_min_pct = 85.0),
        make_screened(35.0, 60.0; coverage_min_pct = 99.0, coverage_max_pct = 99.0)
    )
    for estimate in estimates
        candidates = collect(estimate)
        @test occursin("$(length(estimate)) candidates", sprint(summary, estimate))
        @test occursin(string(estimate.status), sprint(show, estimate))
        @test LineCableModels.TextDisplay.name(typeof(estimate)) == "WireEstimate"
        @test sprint(show, MIME"text/plain"(), estimate; context = :compact=>true) ==
              sprint(show, estimate)
        shown = sprint(show, MIME"text/plain"(), estimate;
            context = IOContext(IOBuffer(), :limit=>false, :displaysize=>(80, 160)))
        @test occursin("target", shown)
        @test occursin("candidates", shown)
        @test occursin(string(estimate.status), shown)
        for reason in estimate.reasons
            @test occursin(reason, shown)
        end
        # Candidates share their display method; the ranked endpoints suffice.
        for candidate in (first(candidates), last(candidates))
            @test occursin("$(candidate.wires) wires", sprint(summary, candidate))
            compact = sprint(show, candidate)
            @test occursin("wires=$(candidate.wires)", compact)
            @test occursin("d=", compact)
            @test sprint(show, MIME"text/plain"(), candidate) == compact
            if hasproperty(candidate, :coverage_pct)
                @test occursin("coverage=", compact)
                @test occursin("%", compact)
            end
        end
        @test collect(estimate) == candidates
    end
end

@testitem "TextDisplay / lazy parametric and uncertainty intent" tags=[:unit] setup=[
    TestFixtures,
] begin
    using Measurements
    using Distributions
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ
    relative=Grid((10.0, 100.0), (1.0, 2.0))
    absolute=Grid((10.0, 100.0), AbsoluteError((0.1, 0.2)))
    calls=Ref(0)
    design=TestFixtures.mv_cable_design()
    space=Gridspace{CableConstantsProblem}(
        temperature->begin
            calls[]+=1
            CableConstantsProblem(design; temperature)
        end,
        (Grid((20.0, 40.0)),)
    )
    nested=Gridspace{CableConstantsProblem}(identity, (space,))
    parametric=ParametricProblem(nested)
    inner=CableConstantsFormulation()
    combinatorial=Combinatorial(inner)
    linear=LinearError(inner; options = (retain_details = true,))
    automatic=MonteCarlo(inner)
    explicit=MonteCarlo(inner; trials = 12, seed = 17, distribution = Uniform(-1, 1),
        return_samples = true, return_histograms = true)
    earth_definition=formula(:default; order = :before, hooks = (contribution = identity,))
    soil_definition=formula(:default; parameters = (tolerance = 1e-6,))
    constants=CableConstants(1e-4, 2e-7, 3e-10)
    results=ParametricResult(combinatorial, [constants])
    empty_results=ParametricResult(combinatorial, typeof(constants)[])
    linear_results=LinearErrorResult(linear, [constants])
    empty_linear=LinearErrorResult(linear, typeof(constants)[])
    mc_results=TestFixtures.cable_monte_carlo_result()
    histogram=UQ.HistogramDensity([1.0, 3.0, 5.0], [0.25, 0.25])
    objects=(
        relative, absolute, AbsoluteError((0.1, 0.2)), PB.UncertainValue(10.0, 0.1),
        Grid((:default, :Ametani2004)), nested, parametric, combinatorial,
        linear, automatic, explicit, results, empty_results,
        linear_results, empty_linear, mc_results, histogram,
        earth_definition, soil_definition, Grid((earth_definition, soil_definition))
    )
    for object in objects
        compact=sprint(show, object)
        @test !occursin('\n', compact)
        @test !isempty(sprint(summary, object))
        @test sprint(show, MIME"text/plain"(), object; context = :compact=>true) == compact
        wide=sprint(show, MIME"text/plain"(), object;
            context = IOContext(IOBuffer(), :limit=>true, :displaysize=>(40, 140)))
        @test !endswith(wide, '\n')
        @test !isempty(wide)
    end
    @test calls[] == 0
    @test occursin("relative error", sprint(show, MIME"text/plain"(), relative))
    @test occursin("%", sprint(show, relative))
    @test occursin("absolute error", sprint(show, MIME"text/plain"(), absolute))
    @test !occursin("%", sprint(show, absolute))
    @test occursin("2 points", sprint(show, nested))
    @test occursin("CableConstantsProblem", sprint(show, nested))
    @test occursin("retained", sprint(show, MIME"text/plain"(), linear))
    @test occursin("DKW-sized", sprint(show, automatic))
    shown=sprint(show, MIME"text/plain"(), explicit)
    @test occursin("12", shown)
    @test occursin("17", shown)
    @test occursin("Uniform", shown)
    @test occursin("samples", shown)
    @test occursin("histograms", shown)
    @test occursin("none", sprint(show, MIME"text/plain"(), empty_results))
    @test occursin("none", sprint(show, MIME"text/plain"(), empty_linear))
    @test occursin("CableConstants", sprint(show, MIME"text/plain"(), linear_results))
    @test occursin("4 per point", sprint(show, MIME"text/plain"(), mc_results))
    @test occursin("2 bins", sprint(show, histogram))
    @test occursin("normalized", sprint(show, MIME"text/plain"(), histogram))
    @test occursin("default", sprint(show, earth_definition))
    @test occursin("before", sprint(show, earth_definition))
    @test occursin("contribution", sprint(show, earth_definition))
    @test occursin("tolerance", sprint(show, soil_definition))
    @test !occursin("{", sprint(show, Grid((earth_definition, soil_definition))))
end

@testitem "TextDisplay / report definitions do not publish or write" tags=[:unit] begin
    using DataFrames
    const RB = LineCableModels.ReportBuilder
    calls = Ref(0)
    illustration = published -> begin
        calls[] += 1
        :illustration
    end
    plain = TableReportDefinition((R, L))
    illustrated = TableReportDefinition((R, L); illustration, clip = false)
    line = RB.LineParametersTableDefinition(
        (@observe(R[:, :, :]),), :base, :kilo, nothing, true)
    mc = RB.MonteCarloTableDefinition(:kilo, nothing, false)
    table = DataFrame(R = [1.0, 2.0], L = [3.0, 4.0])
    original = copy(table)
    artifact = ReportArtifact(table, nothing, nothing)
    illustrated_artifact = ReportArtifact(table, :illustration, :destination)

    mktempdir() do directory
        destination = joinpath(directory, "not-written.xlsx")
        workbook = RB.XLSXWorkbook(destination,
            [
                RB.XLSXSheet("Z(1,1)", ["frequency" "R"; "50" "1"]),
                RB.XLSXSheet("Y(1,1)", ["frequency" "G"; "50" "2"])
            ])
        objects = (
            plain, illustrated, RB.CableConstantsTableDefinition(false), line,
            RB.BenchmarkTableDefinition(false), mc, XLSXReportDefinition(),
            XLSXReportDefinition(file_name = destination, clip = false),
            artifact, illustrated_artifact, workbook
        )
        for object in objects
            compact = sprint(show, object)
            @test !occursin('\n', compact)
            @test !isempty(sprint(summary, object))
            @test sprint(show, MIME"text/plain"(), object; context = :compact=>true) ==
                  compact
            shown = sprint(show, MIME"text/plain"(), object)
            @test !endswith(shown, '\n')
        end
        @test !ispath(destination)
        @test isempty(readdir(directory))
        @test occursin("Z(1,1)", sprint(show, MIME"text/plain"(), workbook))
        @test occursin("Y(1,1)", sprint(show, MIME"text/plain"(), workbook))
        @test occursin("2×2 cells", sprint(show, MIME"text/plain"(), workbook))
    end
    @test calls[] == 0
    @test table == original
    @test occursin("2 observations", sprint(summary, plain))
    @test !occursin("illustration", sprint(show, MIME"text/plain"(), plain))
    @test occursin("illustration", sprint(show, MIME"text/plain"(), illustrated))
    @test occursin("false", sprint(show, MIME"text/plain"(), illustrated))
    @test occursin("kilo", sprint(show, line))
    @test occursin("false", sprint(show, mc))
    @test occursin("default path", sprint(show, XLSXReportDefinition()))
    @test occursin("table=2×2", sprint(show, artifact))
    @test occursin("illustration=none", sprint(show, artifact))
    @test occursin("output=none", sprint(show, artifact))
    @test occursin("illustration=present", sprint(show, illustrated_artifact))
    @test occursin("output=present", sprint(show, illustrated_artifact))
end
