@testitem "Engine / shunt model / strict selection and completed coefficients" tags=[:integration, :importexport, :slow] setup=[SavedResults] begin
    using Logging
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    design = internal_shunt_test_design(count = 4)
    problem = CableConstantsProblem(design)
    default = CableConstantsFormulation()
    equivalent = CableConstantsFormulation(shunt_model = :equivalent)
    nominal = @inferred compute(problem, default)
    @test nominal == compute(problem, equivalent)
    @test details(nominal).data.shunt_model.requested === :equivalent
    @test details(nominal).data.shunt_model.solves == 0
    @test isempty(details(nominal).data.shunt_model.diagnostics)
    @test only(details(nominal).data.shunt_model.domains).terminals == 1:3
    @test occursin("no boundary extraction",
        only(details(nominal).data.shunt_model.domains).message)
    @test isempty(only(E.internal_potential_coefficients([design],
        [LineCableModels.DataModel.flatten(design)], [default])).blocks)
    @test E.formula_id(default.formulas.shunt_model) === :equivalent
    invalid = (formula(:no_such_model), formula(:equivalent; options = (audit = true,)),
        formula(:boundary; options = (audit = true,)),
        formula(:boundary; parameters = (fallback = :silent,)),
        formula(:boundary; parameters = (fallback = :error,)),
        formula(:boundary; parameters = (fallback = :equivalent,)),
        formula(:boundary; options = (integration = (rtol = -1.0,),)),
        formula(:boundary; options = (resolution = (wire = 0,),)),
        formula(:boundary; options = (resolution = (order = 8, quadrature = 8),)))
    for value in invalid
        @test_throws ArgumentError Formulation(shunt_model = value)
        @test_throws ArgumentError CableConstantsFormulation(shunt_model = value)
    end
    # A law that the boundary shunt model cannot represent fails before any solve.
    strict = CableConstantsFormulation(shunt_model = :boundary, insulation_admittance = :lossy)
    @test_throws BoundarySolveError compute(problem, strict)
    # UQ accepts a fixed geometric boundary model. Its failures propagate.
    @test MonteCarlo(strict; trials = 2).inner === strict
    @test LinearError(strict).inner === strict
    @test_throws BoundarySolveError compute(problem, strict)
    quadrature = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),
            integration = (rtol = 1e-15, atol = 0.0, maxevals = 1))))
    logger = Test.TestLogger(; min_level = Logging.Warn)
    approximate = with_logger(logger) do
        compute(problem, quadrature)
    end
    warnings = filter(log -> occursin("logarithmic quadrature target", log.message), logger.logs)
    # One solve reports its unmet QuadGK targets once.
    @test length(warnings) == 1
    @test only(warnings).kwargs[:integrals] > 0 && only(warnings).kwargs[:maxevals] == 1
    @test all(isfinite, approximate.C)
    @test details(approximate).data.shunt_model.requested === :boundary
    @test details(approximate).data.shunt_model.solves == 1
    repeated_logger = Test.TestLogger(; min_level = Logging.Warn)
    unchanged = with_logger(repeated_logger) do
        compute(problem, quadrature)
    end
    @test unchanged.C == approximate.C
    @test details(unchanged).data.shunt_model.requested === :boundary
    @test all(log -> !occursin("replaced by the equivalent annular layer", log.message), repeated_logger.logs)
    lossy = CableConstantsFormulation(shunt_model = :boundary, insulation_admittance = :lossy)
    @test_throws BoundarySolveError compute(problem, lossy)
    boundary = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),)))
    blueprint = @inferred LineCableModels.DataModel.flatten(design)
    shunt = only(E.internal_potential_coefficients([design], [blueprint], [boundary]))
    @test length(shunt.blocks) == 1
    @test shunt.details.solves == 1
    @test only(E.internal_potential_coefficients([design],
        [LineCableModels.DataModel.flatten(design, Float32)], [boundary])).blocks isa
          Vector{E.InternalCapacitance{Float32}}
    first_result = compute(problem, boundary)
    explicit_lossless = CableConstantsFormulation(
        shunt_model = boundary.formulas.shunt_model,
        insulation_admittance = :lossless, semicon_admittance = :lossless)
    @test compute(problem, explicit_lossless) == first_result
    workspace = E.CableConstantsWorkspace(problem, boundary, E.LocalCableData(blueprint), shunt)
    @test workspace.plan.shunt.blocks[1].C === shunt.blocks[1].C
    @test E._solve!(workspace, problem, boundary) == first_result
    @test details(first_result).data.shunt_model.requested === :boundary
    @test keys(only(details(first_result).data.shunt_model.diagnostics)) ==
          (:log_moment_error, :condition, :sectors, :unknowns, :equations, :level)
    @test typeof(first_result) === typeof(nominal)
    @test compute(problem, boundary) == first_result
    @test compute(CableConstantsProblem(design; frequency = 60), boundary).C ==
          first_result.C
    changed = CableConstantsProblem(internal_shunt_test_design(count = 4, epsilon = 4.0))
    @test compute(changed, boundary).C != first_result.C
    @test compute(problem, default) == nominal
    for id in E.ShuntModel.formulas(E.ShuntModel.Formula)
        selected = E.ShuntModel.Formula(id)
        @test occursin("id=" * string(formula_id(selected)), sprint(show, selected))
        @test occursin(string(formula_id(selected)), sprint(show, MIME"text/plain"(), selected))
    end
    @test E.BoundarySolveError === E.ShuntModel.BoundarySolveError
    @test which(show, (IO, MIME"text/plain", BoundarySolveError)).module === E.ShuntModel
    for formulation in
        (boundary, Formulation(shunt_model = boundary.formulas.shunt_model))
        @test occursin("boundary", description(formulation, C))
    end
    for source in (first_result, nominal)
        @test SavedResults.same(SavedResults.round_trip(source), source)
    end
end
