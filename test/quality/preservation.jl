# Preservation locks. `preservation.toml` records counts per key, in tables that
# declare their direction: `[inferred]` is a floor and `[jet]` a ceiling. A check fails
# when a live count passes its limit, and also when the live count improves on the
# table, so the table always records the current state. `test/tools/baseline_ratchet.jl
# REF` checks that floors only rise and ceilings only fall across commits.
@testmodule PreservationLocks begin
    using TestItemRunner, TOML
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "runner.jl"))
    const T = ValidationTestRunner
    const REPOSITORY = dirname(dirname(@__DIR__))
    const TABLES = joinpath(@__DIR__, "preservation.toml")

    table(name) = Dict{String, Int}(get(TOML.parsefile(TABLES), name, Dict{String, Any}()))
    isfloor(name) = TOML.parsefile(TABLES)["directions"][name] == "floor"

    # Keys whose live count passes the limit (`broken`) and keys whose live count
    # improves on it (`stale`). A missing key counts as zero.
    function compare(live, recorded; floor::Bool)
        broken, stale = String[], String[]
        for key in sort!(collect(union(keys(live), keys(recorded))))
            now, limit = get(live, key, 0), get(recorded, key, 0)
            now == limit && continue
            push!((now < limit) == floor ? broken : stale, "$key: $now (recorded $limit)")
        end
        return (; broken, stale)
    end

    report(lines) = (foreach(println, lines); lines)

    # L1. The `@inferred` uses in a Julia source, qualified or not.
    function inferred_uses(text, filename = "none")
        uses = Ref(0)
        visit(x) = if x isa Expr
            if Meta.isexpr(x, :macrocall)
                name = x.args[1]
                name isa Expr && Meta.isexpr(name, :.) && (name = name.args[end])
                name isa QuoteNode && (name = name.value)
                name === Symbol("@inferred") && (uses[] += 1)
            end
            foreach(visit, x.args)
        end
        visit(Meta.parseall(text; filename))
        return uses[]
    end

    # The test files with at least one `@inferred` use, and their counts.
    function inferred_counts(root)
        counts = Dict{String, Int}()
        for file in TestItemRunner.find_test_files(root)
            startswith(file, joinpath(root, "test", "")) || continue
            uses = inferred_uses(read(file, String), file)
            uses > 0 && (counts[T.relative(file, root)] = uses)
        end
        return counts
    end
end

@testitem "Quality / preservation / L1 @inferred floor" tags=[:quality] setup=[PreservationLocks] begin
    P = PreservationLocks
    result = P.compare(P.inferred_counts(P.REPOSITORY), P.table("inferred");
        floor = P.isfloor("inferred"))
    # A file below its floor lost `@inferred` assertions.
    @test P.report(result.broken) == String[]
    # A file above its floor raises the floor in the same change.
    @test P.report(result.stale) == String[]
end

@testitem "Quality / preservation / JET ceilings on the frequency-loop kernels" tags=[:quality] setup=[PreservationLocks] begin
    using JET
    P = PreservationLocks
    E, M, C = LineCableModels.Engine, LineCableModels.ModalAnalysis, LineCableModels.Commons
    include(joinpath(P.REPOSITORY, "test", "support", "scenarios.jl"))
    using .CurrentScenarios: line_parameters_problem, three_phase_system

    # The default coaxial calculation on a cable system, run once so that every buffer
    # holds the state of the frequency loop. Each kernel is analysed with the concrete
    # arguments of its call in `_solve!`.
    problem = line_parameters_problem(three_phase_system(); frequencies = [50.0, 1000.0])
    formulation = Formulation()
    blueprint = only(E.flatten(LineCableModelsCoaxial(), problem.system.designs,
        eltype(problem), [formulation]))
    execution = E.computation_options(LineCableModelsCoaxial, ComputationOptions())
    workspace = E.LineParametersWorkspace(problem, formulation, execution,
        E.lineinput(problem, blueprint))
    E._solve!(workspace, formulation)
    buffers, invariants, input = workspace.buffers, workspace.invariants, workspace.input
    jω = input.jω[1]
    found = Dict{String, Vector{Any}}()
    found["coaxial _solve!"] = JET.get_reports(
        @report_opt target_modules=(LineCableModels,) E._solve!(workspace, formulation))
    found["cable_impedance!"] = JET.get_reports(
        @report_opt target_modules=(LineCableModels,) E.cable_impedance!(buffers.Zprimitive,
            input.cable, buffers.rho_cond, formulation.methods, jω; workspace))
    found["cable_potential!"] = JET.get_reports(
        @report_opt target_modules=(LineCableModels,) E.cable_potential!(buffers.Pprimitive,
            input.cable, buffers.dielectric_admittivity, jω, buffers.layer_coefficients,
            buffers.coefficients, buffers.tails))
    found["earth!"] = JET.get_reports(
        @report_opt target_modules=(LineCableModels,) E.earth!(workspace, 1,
            invariants.earth_calculations, buffers.earth_materials))
    found["reduce_line_matrices!"] = JET.get_reports(
        @report_opt target_modules=(LineCableModels,) C.reduce_line_matrices!(
            view(buffers.Zout, :, :, 1), view(buffers.Yout, :, :, 1), buffers.Zprimitive,
            buffers.Pprimitive, jω, invariants.plan, buffers.reduction))
    phase = compute(problem, formulation)
    for id in (:chrysochos2014, :vieira2026, :wedepohl1996)
        selected = ModalAnalysisFormulation(id).formula
        modal = M.ModalAnalysisWorkspace(phase, selected)
        found["decompose! $id"] = JET.get_reports(
            @report_opt target_modules=(LineCableModels,) M.decompose!(
                M.allocation_selector(selected), modal, M.formula_parameters(selected),
                M.formulation_options(selected)))
    end
    for (kernel, reports) in sort!(collect(found); by = first), report in reports
        println(kernel, ": ", first(split(sprint(show, report), '\n')))
    end
    result = P.compare(Dict(k => length(v) for (k, v) in found), P.table("jet");
        floor = P.isfloor("jet"))
    # A kernel above its ceiling has a new optimization report.
    @test P.report(result.broken) == String[]
    # A kernel below its ceiling lowers the ceiling in the same change.
    @test P.report(result.stale) == String[]
end

@testitem "Quality / preservation / negative controls" tags=[:quality] setup=[PreservationLocks] begin
    P = PreservationLocks
    # `@inferred` counts bare and qualified uses, nested ones too, but not text.
    source = """
        @inferred f(1)
        @test (@inferred g()) == 1
        Test.@inferred h(2)
        # @inferred in a comment
        text = "@inferred in a string"
        """
    @test P.inferred_uses(source) == 3
    @test P.inferred_uses("x = 1") == 0

    # A floor breaks below its limit and a ceiling above it; improvements are stale.
    floor = P.compare(Dict("a" => 2, "b" => 5, "new" => 1), Dict("a" => 3, "b" => 4, "gone" => 1);
        floor = true)
    @test floor.broken == ["a: 2 (recorded 3)", "gone: 0 (recorded 1)"]
    @test floor.stale == ["b: 5 (recorded 4)", "new: 1 (recorded 0)"]
    ceiling = P.compare(Dict("a" => 2, "b" => 5), Dict("a" => 3, "b" => 4); floor = false)
    @test ceiling.broken == ["b: 5 (recorded 4)"]
    @test ceiling.stale == ["a: 2 (recorded 3)"]
    @test P.compare(Dict("a" => 1), Dict("a" => 1); floor = true) == (; broken = String[], stale = String[])

    # The ratchet follows each declared direction across commits.
    ratchet = Module(:BaselineRatchet)
    Base.include(ratchet, joinpath(P.REPOSITORY, "test", "tools", "baseline_ratchet.jl"))
    mktempdir() do repository
        git(arguments...) = run(pipeline(Base.invokelatest(ratchet.git, repository,
            "-c", "user.name=Probe", "-c", "user.email=probe@example.invalid",
            arguments...); stdout = devnull, stderr = devnull))
        function locks(inferred, jet; directions = "inferred = \"floor\"\njet = \"ceiling\"\n")
            path = joinpath(repository, ratchet.PRESERVATION)
            mkpath(dirname(path))
            write(path, "[directions]\n" * directions *
                "[inferred]\n" * join(("\"$k\" = $v\n" for (k, v) in inferred)) *
                "[jet]\n" * join(("\"$k\" = $v\n" for (k, v) in jet)))
        end
        moved() = Base.invokelatest(ratchet.moved, repository, "HEAD", ratchet.PRESERVATION).lines
        git("init", "-q")
        mkpath(joinpath(repository, "test"))
        write(joinpath(repository, "test", "a.jl"), "@inferred f()\n")
        locks(["test/a.jl" => 3], ["kernel" => 2])
        git("add", "-A")
        git("commit", "-q", "--no-verify", "--no-gpg-sign", "-m", "start")
        locks(["test/a.jl" => 4, "test/new.jl" => 1], ["kernel" => 1])
        @test moved() == String[]
        locks(["test/a.jl" => 2], ["kernel" => 2])
        @test moved() == ["inferred | test/a.jl: 2 (3 at HEAD)"]
        locks(Pair{String, Int}[], ["kernel" => 3, "other" => 1])
        @test moved() == ["inferred | test/a.jl: removed (3 at HEAD)",
            "jet | kernel: 3 (2 at HEAD)", "jet | other: 1 (absent at HEAD)"]
        # A floor key follows a file renamed by git.
        git("mv", "test/a.jl", "test/b.jl")
        locks(["test/b.jl" => 3], ["kernel" => 2])
        @test moved() == String[]
        locks(["test/b.jl" => 3], ["kernel" => 2]; directions = "inferred = \"floor\"\n")
        @test_throws ErrorException moved()
    end
end
