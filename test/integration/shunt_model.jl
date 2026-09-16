@testitem "Engine / shunt model / selection, fallback and blueprint coefficients" tags=[:integration] begin
    E = LineCableModels.Engine
    IE = LineCableModels.ImportExport
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    design = internal_shunt_test_design(count = 4)
    problem = CableConstantsProblem(design)
    default = CableConstantsFormulation()
    coaxial = CableConstantsFormulation(shunt_model = :coaxial)
    nominal = @inferred compute(problem, default)
    @test nominal == compute(problem, coaxial)
    @test details(nominal).data.shunt_model.effective === :coaxial
    @test details(nominal).data.shunt_model.solves == 0
    @test isempty(E.flatten(LineCableModelsCoaxial(), design, default).shunt)
    @test E.formula_id(default.methods.shunt_model) === :coaxial
    invalid = (formula(:no_such_model), formula(:coaxial; options = (audit = true,)),
        formula(:boundary; parameters = (fallback = :silent,)),
        formula(:boundary; options = (integration = (rtol = -1.0,),)),
        formula(:boundary; options = (resolution = (wire = 0,),)),
        formula(:boundary; options = (resolution = (order = 8, quadrature = 8),)))
    for value in invalid
        @test_throws ArgumentError Formulation(shunt_model = value)
        @test_throws ArgumentError CableConstantsFormulation(shunt_model = value)
    end
    # Exhaust a declared numerical budget before allocating a dense matrix.
    budget = (resolution = (wire = 100_000,),)
    strict = CableConstantsFormulation(shunt_model = formula(:boundary; options = budget))
    @test_throws BoundarySolveError compute(problem, strict)
    fallback = CableConstantsFormulation(shunt_model = formula(:boundary;
        parameters = (fallback = :coaxial,), options = budget))
    result = @test_logs (:warn, r"replaced by coaxial") compute(problem, fallback)
    report = details(result).data.shunt_model
    @test result == nominal
    @test report.requested === :boundary && report.effective === :coaxial
    @test only(report.domains).reason === :budget
    @test occursin("budget", only(report.domains).message)
    @test_throws ArgumentError MonteCarlo(fallback; trials = 2)
    @test_throws ArgumentError LinearError(fallback)
    quadrature = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),
            integration = (rtol = 1e-15, atol = 0.0, maxevals = 1))))
    error = try
        compute(problem, quadrature)
    catch exception
        exception
    end
    @test error isa BoundarySolveError
    @test error.category === :quadrature
    @test error.context.stage === :assembly
    @test error.context.design == 1
    @test error.context.evaluations > 0
    @test error.context.estimate > error.context.tolerance
    lossy = CableConstantsFormulation(shunt_model = :boundary, insulation_admittance = :lossy)
    @test_throws BoundarySolveError compute(problem, lossy)
    boundary = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),)))
    blueprint = @inferred E.flatten(LineCableModelsCoaxial(), design, boundary)
    @test length(blueprint.shunt) == 1
    @test blueprint.shunt_details.solves == 1
    @test E.flatten(LineCableModelsCoaxial(), design, Float32, boundary) isa
          E.CableBlueprint{Float32}
    first_result = compute(problem, boundary)
    explicit_lossless = CableConstantsFormulation(
        shunt_model = boundary.definitions.shunt_model,
        insulation_admittance = :lossless, semicon_admittance = :lossless)
    @test compute(problem, explicit_lossless) == first_result
    workspace = E.CableConstantsWorkspace(problem, boundary, blueprint)
    @test workspace.cable.shunt[1].C === blueprint.shunt[1].C
    @test E._solve!(workspace, problem, boundary) == first_result
    @test details(first_result).data.shunt_model.effective === :boundary
    @test only(details(first_result).data.shunt_model.diagnostics).boundary_residual === nothing
    @test typeof(first_result) === typeof(nominal)
    @test compute(problem, boundary) == first_result
    @test compute(CableConstantsProblem(design; frequency = 60), boundary).C ==
          first_result.C
    changed = CableConstantsProblem(internal_shunt_test_design(count = 4, epsilon = 4.0))
    @test compute(changed, boundary).C != first_result.C
    @test compute(problem, default) == nominal
    for id in E.ShuntModel.formulas()
        selected = E.ShuntModel.Formula(id)
        @test occursin("id=" * string(formula_id(selected)), sprint(show, selected))
        @test occursin(string(formula_id(selected)), sprint(show, MIME"text/plain"(), selected))
    end
    @test which(show, (IO, MIME"text/plain", BoundarySolveError)).module === E
    @test occursin("quadrature", sprint(show, MIME"text/plain"(), error))
    for formulation in
        (boundary, Formulation(shunt_model = boundary.definitions.shunt_model))
        record=IE.deserialize_value(Val(:formulation), NamedTuple(formulation))
        @test formula_id(record, C) == formula_id(formulation, C)
        @test occursin("boundary", description(formulation, C))
    end
    for source in (first_result, result)
        restored=IE.deserialize_value(IE.serialize_value(source))
        @test restored == source
        @test details(restored).data.shunt_model.effective ===
              details(source).data.shunt_model.effective
        @test only(details(restored).data.shunt_model.domains).reason ===
              only(details(source).data.shunt_model.domains).reason
    end
end
