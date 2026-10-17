@testitem "LineModels / completed segment and fitting selection composition" tags=[:unit, :parametric] setup=[TestFixtures] begin
    import LineCableModels.LineModels as LM

    # Record each scalar fallback so the finite traversal can complete without the extension.
    struct FittingOutcome{F <: LineModelFormulation} <: AbstractFormulation
        selected::F
    end
    struct MissingFittingExtension{S} <: AbstractCoreResult
        segment::S
        identifier::Symbol
        options::NamedTuple
        message::String
    end
    function LineCableModels.compute(problem::LineModelProblem, formulation::FittingOutcome;
            options::Union{NamedTuple, ComputationOptions}=ComputationOptions())
        outcome = try
            selected = formulation.selected.formula
            invoke(LM.fit, Tuple{LM.Formula{formula_id(selected)}, Vararg{Any}},
                selected, problem, selected.parameters, selected.options)
        catch exception
            exception isa ArgumentError || rethrow()
            exception
        end
        return MissingFittingExtension(problem.segment, formula_id(formulation.selected.formula),
            formulation_options(formulation.selected.formula).data, sprint(showerror, outcome))
    end

    phase = TestFixtures.two_conductor_results()
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation())
    segments = ParametricResult(Combinatorial(ModalAnalysisFormulation()),
        collect(PropagationParameters(modal; line_length=Grid((500.0, 1000.0)))))
    problems = Gridspace{LineModelProblem}(segments)
    formulations = LineModelFormulation(Grid([:ulm, :fdcm]);
        options=(yc=(npoles=8,), h=(npoles=8,), passivity=false))
    @test length(problems) == length(segments) == 2
    @test [point.segment for point in problems] == collect(segments)
    @test formulations isa Gridspace{LineModelFormulation}
    @test length(formulations) == 2
    @test [formula_id(point.formula) for point in formulations] == [:ulm, :fdcm]

    outcomes = Gridspace{FittingOutcome}(FittingOutcome, (formulations,))
    run = compute(problems, outcomes)
    @test run isa ParametricResult
    @test length(run) == 4
    @test run.axes.problems === problems
    @test [value.selected for value in run.axes.formulations] == collect(formulations)
    for j in 1:length(formulations), i in 1:length(segments)
        point = run[i, j]
        @test point === run[i + (j - 1) * length(segments)]
        @test point.segment === segments[i]
        @test point.options == formulation_options(run.axes.formulations[j].selected.formula).data
        @test point.identifier == formula_id(run.axes.formulations[j].selected.formula)
        @test point.message == "ArgumentError: load UniversalLineModel to fit :$(point.identifier)"
    end
end

@testitem "LineModels / completion details across option values" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels.LineModels as LM

    struct FixedFit{O <: FormulationOptions} <: AbstractFormulation
        parameters::NamedTuple
        options::O
    end
    LM.formula_id(::FixedFit) = :fixed_fit
    LM.initialize_buffers(::FixedFit, ::Type, input, plan, buffers) = buffers
    function LM.fit(::FixedFit, workspace, parameters, options)
        n = length(workspace.input.coordinates)
        admittance = (poles=ComplexF64[-1], residues=ones(ComplexF64, n, n, 1), direct=zeros(n, n))
        groups = [(delay=1e-6, poles=ComplexF64[-2], residues=ones(ComplexF64, n, n, 1))]
        errors = (rms=zeros(n, n), max=zeros(n, n), f_max=fill(first(workspace.input.frequencies), n, n))
        passive = (checked=false, passive=nothing, bands=Tuple{Float64, Float64}[],
            worst=nothing, f_worst=nothing, perturbation=nothing)
        return (WidebandModel(admittance, groups), (admittance=errors, propagation=errors),
            passive, (targets=nothing, achieved=nothing), (;))
    end

    phase = TestFixtures.two_conductor_results()
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation())
    problem = LineModelProblem(PropagationParameters(modal; line_length=1000.0))
    controls = ((yc=(npoles=8,), h=(npoles=8,), passivity=false),
        (yc=(npoles=10,), h=(npoles=12,), passivity=(;)))
    models = map(controls) do options
        compute(problem, LineModelFormulation(FixedFit((;), FormulationOptions(options))))
    end
    @test all(model -> model isa LineModel, models)
    @test isconcretetype(typeof(first(models)))
    @test typeof(first(models)) == typeof(last(models))
    @test details(first(models)).data.options == first(controls)
    @test details(last(models)).data.options == last(controls)
end
