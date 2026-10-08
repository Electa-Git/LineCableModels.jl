# The preservation locks record counts per key in `preservation.toml`. Each table
# declares its direction. `[inferred]` is a floor, and `[jet]` and `[allocations]` are
# ceilings. A check fails when a live count passes its limit. It also fails when the
# live count improves on the table, so that the table records the current state.
# Allocation counts may differ from the table by their allowance, and recorded bytes fail
# only on increases beyond their margin. The measured ceilings hold for the Julia
# version and the `Manifest.toml` recorded in `[environment]`. Across commits,
# `test/tools/baseline_ratchet.jl` checks that floors only rise and ceilings only fall.
@testmodule PreservationLocks begin
    using TestItemRunner, TOML, SHA
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

    # The environment of the measured tables: the Julia version and the committed
    # `Manifest.toml`, which the CI quality job instantiates.
    const MANIFEST = joinpath(REPOSITORY, "Manifest.toml")
    environment() = Dict("julia" => string(VERSION),
        "manifest" => isfile(MANIFEST) ? bytes2hex(sha256(read(MANIFEST))) : "absent")
    describe(e) = "Julia $(e["julia"]) and Manifest $(first(e["manifest"], 12))"

    # The measured tables hold for the recorded environment only.
    function check_environment(recorded = TOML.parsefile(TABLES)["environment"],
            running = environment())
        recorded == running || error("recorded with $(describe(recorded)), running " *
            "$(describe(running)): re-record [jet] and [allocations]")
        return nothing
    end

    # Allocation ceilings. Measurements stores a partial derivative only when it is
    # nonzero, but whether a round-off derivative is exactly zero depends on the machine's
    # last bits. So a count may differ from the table by the scenario's allowance, the number of
    # uncertain real scalars that its result publishes. The allocation tool computes it
    # live. A plain-number scenario has allowance zero and compares exactly. Bytes fail
    # above the recorded minimum plus `HEADROOM` for the runtime's byte accounting and the
    # allowance times the bytes of one derivative entry.
    const HEADROOM = 512
    function compare_allocations(live, recorded, allowances, derivative)
        broken, stale = String[], String[]
        for key in sort!(collect(union(keys(live), keys(recorded))))
            now, limit = get(live, key, 0), get(recorded, key, 0)
            allowance = get(allowances, first(rsplit(key, " | "; limit = 2)), 0)
            if endswith(key, "| allocations")
                line = "$key: $now (recorded $limit ± $allowance)"
                now > limit + allowance && push!(broken, line)
                now < limit - allowance && push!(stale, line)
            elseif endswith(key, "| bytes")
                margin = HEADROOM + allowance*derivative
                now > limit + margin && push!(broken, "$key: $now (recorded $limit + $margin)")
            end
        end
        return (; broken, stale)
    end

    const ALLOCATIONS = joinpath(REPOSITORY, "test", "tools", "allocations.jl")

    # The `@inferred` floor: the `@inferred` uses in a Julia source, qualified or not.
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

@testitem "Quality / preservation / @inferred floor" tags=[:quality] setup=[PreservationLocks] begin
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
    P.check_environment()
    E, M, C = LineCableModels.Engine, LineCableModels.ModalAnalysis, LineCableModels.Commons
    include(joinpath(P.REPOSITORY, "test", "support", "scenarios.jl"))
    using .CurrentScenarios: line_parameters_problem, three_phase_system

    # The default coaxial computation on a cable system, run once so that every buffer
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
    buffers, plan, input = workspace.buffers, workspace.plan, workspace.input
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
            plan.earth.calculations, buffers.earth.calculations))
    found["reduce_line_matrices!"] = JET.get_reports(
        @report_opt target_modules=(LineCableModels,) C.reduce_line_matrices!(
            view(buffers.Zout, :, :, 1), view(buffers.Yout, :, :, 1), buffers.Zprimitive,
            buffers.Pprimitive, jω, plan.reduction, buffers.reduction))
    phase = compute(problem, formulation)
    for id in (:chrysochos2014, :vieira2026, :wedepohl1996)
        selected = ModalAnalysisFormulation(id).formula
        modal = M.ModalAnalysisWorkspace(phase, selected)
        found["decompose! $id"] = JET.get_reports(
            @report_opt target_modules=(LineCableModels,) M.decompose!(
                selected, modal, M.formula_parameters(selected),
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

@testitem "Quality / preservation / allocation ceilings" tags=[:quality] setup=[PreservationLocks] begin
    using TOML
    P = PreservationLocks
    P.check_environment()
    # A fresh process runs the corpus in a fixed order: the counts depend on what the
    # process computed before.
    output, errors = IOBuffer(), IOBuffer()
    command = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) $(P.ALLOCATIONS)`
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = errors))
    text = String(take!(output))
    print(text)
    # A nonzero exit reports allocation counts that varied across calls.
    success(process) || print(String(take!(errors)))
    @test success(process)
    measured = TOML.parse(text)
    live = Dict{String, Int}(measured["allocations"])
    result = P.compare_allocations(live, P.table("allocations"),
        Dict{String, Int}(measured["allowances"]), measured["derivative"]["bytes"])
    # A row above its ceiling allocates more.
    @test P.report(result.broken) == String[]
    # A row whose counts fell lowers its ceiling, and re-records its bytes, in the same change.
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

    # A floor breaks below its limit and a ceiling above it. Improvements are stale.
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
        # An `@inferred` floor can fall when its lines move to another file, unchanged
        # apart from indentation. A deleted line, an edited line, a line re-indented in the
        # same file and a drop beyond the moved lines fail.
        source, target = joinpath(repository, "test", "a.jl"), joinpath(repository, "test", "c.jl")
        write(source, "")
        write(target, "    @inferred f()\n")
        locks(["test/a.jl" => 2, "test/c.jl" => 1], ["kernel" => 2])
        @test moved() == String[]
        locks(["test/c.jl" => 1], ["kernel" => 2])
        @test moved() == ["inferred | test/a.jl: removed (3 at HEAD)"]
        locks(["test/a.jl" => 2, "test/c.jl" => 1], ["kernel" => 2])
        write(target, "@inferred f(1)\n")
        @test moved() == ["inferred | test/a.jl: 2 (3 at HEAD)"]
        rm(target)
        @test moved() == ["inferred | test/a.jl: 2 (3 at HEAD)"]
        write(source, "    @inferred f()\n")
        @test moved() == ["inferred | test/a.jl: 2 (3 at HEAD)"]
        write(source, "@inferred f()\n")
        # A floor key follows a file renamed by git.
        git("mv", "test/a.jl", "test/b.jl")
        locks(["test/b.jl" => 3], ["kernel" => 2])
        @test moved() == String[]
        locks(["test/b.jl" => 3], ["kernel" => 2]; directions = "inferred = \"floor\"\n")
        @test_throws ErrorException moved()
        # A measured ceiling may rise only together with its recorded environment.
        function measured(jet, allocations, julia; manifest = "a")
            write(joinpath(repository, ratchet.PRESERVATION),
                "[environment]\njulia = \"$julia\"\nmanifest = \"$manifest\"\n" *
                "[directions]\ninferred = \"floor\"\n" *
                "jet = \"ceiling\"\nallocations = \"ceiling\"\n[inferred]\n\"test/b.jl\" = 3\n" *
                "[jet]\nkernel = $jet\n[allocations]\n\"s | 2 frequencies | allocations\" = $allocations\n")
        end
        measured(2, 10, "1.12.7")
        git("add", "-A")
        git("commit", "-q", "--no-verify", "--no-gpg-sign", "-m", "measured")
        measured(3, 11, "1.12.7")
        @test moved() == ["allocations | s | 2 frequencies | allocations: 11 (10 at HEAD)",
            "jet | kernel: 3 (2 at HEAD)"]
        measured(3, 11, "1.12.8")
        @test moved() == String[]
        measured(3, 11, "1.12.7"; manifest = "b")
        @test moved() == String[]
        # A floor stays a floor across a version change.
        write(joinpath(repository, ratchet.PRESERVATION),
            replace(read(joinpath(repository, ratchet.PRESERVATION), String),
                "\"test/b.jl\" = 3" => "\"test/b.jl\" = 2"))
        @test moved() == ["inferred | test/b.jl: 2 (3 at HEAD)"]
        # Each allocation key is a corpus scenario. New scenarios pass while every scenario
        # at HEAD keeps its key, and a key that replaces another is a rename that fails.
        function scenarios(rows)
            write(joinpath(repository, ratchet.PRESERVATION),
                "[environment]\njulia = \"1.12.7\"\nmanifest = \"a\"\n" *
                "[directions]\ninferred = \"floor\"\n" *
                "jet = \"ceiling\"\nallocations = \"ceiling\"\n[inferred]\n\"test/b.jl\" = 3\n" *
                "[jet]\nkernel = 2\n[allocations]\n" * join(("\"$k\" = $v\n" for (k, v) in rows)))
        end
        scenarios(["s | 2 frequencies | allocations" => 10, "t | 2 frequencies | allocations" => 7])
        @test moved() == String[]
        scenarios(["t | 2 frequencies | allocations" => 7])
        @test moved() == ["allocations | t | 2 frequencies | allocations: 7 (absent at HEAD)"]
    end

    # The measured tables refuse another Julia version or another Manifest.
    recorded = Dict("julia" => "1.12.7", "manifest" => "0123456789abcdef")
    @test P.check_environment(recorded, copy(recorded)) === nothing
    @test_throws ErrorException("recorded with Julia 1.12.7 and Manifest 0123456789ab, running " *
        "Julia 1.12.8 and Manifest 0123456789ab: re-record [jet] and [allocations]") P.check_environment(
        recorded, merge(recorded, Dict("julia" => "1.12.8")))
    @test_throws ErrorException P.check_environment(recorded, merge(recorded, Dict("manifest" => "fedcba")))

    # A scenario without uncertain results compares exactly in both directions. A scenario
    # with uncertain results can move by its allowance in objects. In bytes, it can exceed
    # the headroom by the allowance times one derivative entry.
    recorded = Dict("s | 2 frequencies | allocations" => 10, "s | 2 frequencies | bytes" => 1000)
    row(n, b) = Dict("s | 2 frequencies | allocations" => n, "s | 2 frequencies | bytes" => b)
    exact, derivative = Dict("s | 2 frequencies" => 0), 48
    compare(n, b, allowances = exact) = P.compare_allocations(row(n, b), recorded, allowances, derivative)
    @test compare(10, 1000 + P.HEADROOM) == (; broken = String[], stale = String[])
    @test compare(10, 900) == (; broken = String[], stale = String[])
    @test compare(10, 1001 + P.HEADROOM).broken ==
        ["s | 2 frequencies | bytes: 1513 (recorded 1000 + 512)"]
    @test compare(11, 1000).broken == ["s | 2 frequencies | allocations: 11 (recorded 10 ± 0)"]
    @test compare(9, 1000).stale == ["s | 2 frequencies | allocations: 9 (recorded 10 ± 0)"]
    # A scenario that the allowances do not list compares exactly.
    @test P.compare_allocations(row(11, 1000), recorded, Dict{String, Int}(), derivative).broken ==
        ["s | 2 frequencies | allocations: 11 (recorded 10 ± 0)"]
    uncertain = Dict("s | 2 frequencies" => 3)
    margin = P.HEADROOM + 3*derivative
    @test compare(13, 1000 + margin, uncertain) == (; broken = String[], stale = String[])
    @test compare(7, 1000, uncertain) == (; broken = String[], stale = String[])
    @test compare(14, 1000, uncertain).broken == ["s | 2 frequencies | allocations: 14 (recorded 10 ± 3)"]
    @test compare(6, 1000, uncertain).stale == ["s | 2 frequencies | allocations: 6 (recorded 10 ± 3)"]
    @test compare(10, 1001 + margin, uncertain).broken ==
        ["s | 2 frequencies | bytes: $(1001 + margin) (recorded 1000 + $margin)"]

    # The allocation tool reports object counts that vary across calls, and records the
    # minimum bytes.
    tool = Module(:AllocationTool)
    Base.include(tool, P.ALLOCATIONS)
    rows(objects, bytes) = Base.invokelatest(tool.scenario_rows, "s", 2, (; objects, bytes))
    @test rows([5, 5, 5], [12, 10, 11]) ==
        (["s | 2 frequencies | allocations" => 5, "s | 2 frequencies | bytes" => 10], nothing)
    @test last(rows([5, 6, 5], [10, 10, 10])) ==
        "s | 2 frequencies: allocations [5, 6, 5] across calls"
    # An allowance counts uncertain real and imaginary parts by type. An exact uncertain
    # value counts, and so does a plain entry that an uncertain array promotes. One
    # derivative entry takes a positive number of bytes.
    in_tool(f, arguments...) = Base.invokelatest(getfield(tool, f), arguments...)
    exact_uncertain = in_tool(:measurement, 1.0, 0.0)
    @test in_tool(:allowance, complex(exact_uncertain, exact_uncertain)) == 2
    @test in_tool(:allowance, [complex(exact_uncertain, exact_uncertain), 1.0 + 2.0im]) == 4
    @test in_tool(:allowance, zeros(ComplexF64, 3)) == 0
    @test in_tool(:derivative_bytes) > 0

    # The equivalence check reads declarations and renames, classifies differences by kind
    # and finds the revision's own files.
    equivalence = Module(:EquivalenceTool)
    Base.include(equivalence, joinpath(P.REPOSITORY, "test", "tools", "equivalence.jl"))
    in_equivalence(f, arguments...; keywords...) =
        Base.invokelatest(getfield(equivalence, f), arguments...; keywords...)
    mktempdir() do directory
        file = joinpath(directory, "allowed")
        write(file, "# Declared differences\nrename | Grammar => Commons\ncoaxial | | value\n" *
            "modal | result.Z | type  # a path prefix\n\n")
        @test in_equivalence(:declarations, file) == (; renames = ["Grammar" => "Commons"], declared = [
            (; scenario = "coaxial", prefix = "", kind = :value),
            (; scenario = "modal", prefix = "result.Z", kind = :type)])
        for malformed in ("coaxial | result", "coaxial | result | values", "rename | Grammar",
                "rename | Grammar => Commons.Types")
            write(file, malformed * "\n")
            @test_throws ErrorException in_equivalence(:declarations, file)
        end
    end
    # A rename changes whole identifiers in type names, path segments and Symbol values.
    nodes = ["result.Grammar.x" => ("Grammar.FormulaDefinition{:Grammar}", ""),
        "result.y" => ("Symbol", ":Grammar"), "result.z" => ("String", "\"Grammar\""),
        "result.w" => ("MyGrammar.GrammarX", ""), "result.eltype" => ("Type", "Grammar.Point")]
    used = falses(2)
    @test in_equivalence(:renamed, nodes, ["Grammar" => "Commons", "Absent" => "Present"], used) == [
        "result.Commons.x" => ("Commons.FormulaDefinition{:Commons}", ""),
        "result.y" => ("Symbol", ":Commons"), "result.z" => ("String", "\"Grammar\""),
        "result.w" => ("MyGrammar.GrammarX", ""), "result.eltype" => ("Type", "Commons.Point")]
    # The second rename is unused, and the report says so.
    @test used == [true, false]
    # Differences come in recorded order with their kind. A declaration covers its kind only.
    before = ["a" => ("Float64", "01"), "b" => ("Int64", "1"), "gone" => ("Int64", "1")]
    after = ["a" => ("Float64", "02"), "b" => ("Int32", "1"), "added" => ("Int64", "1")]
    changes = in_equivalence(:differences, before, after)
    @test first.(changes) == ["a", "b", "added", "gone"]
    @test [in_equivalence(:kind, c) for c in changes] == [:value, :type, :path, :path]
    @test isempty(in_equivalence(:differences, before, before))
    value = (; scenario = "s", prefix = "", kind = :value)
    @test in_equivalence(:covers, value, "s", changes[1])
    @test !in_equivalence(:covers, value, "s", changes[2])
    @test in_equivalence(:covers, (; scenario = "s", prefix = "b", kind = :type), "s", changes[2])
    @test !in_equivalence(:covers, value, "other", changes[1])
    # The timing comparison pins both workers to one performance core of a hybrid machine,
    # never CPU 0 or its sibling. It grades each scenario and sets its exit status.
    timing = Module(:TimingTool)
    Base.include(timing, joinpath(P.REPOSITORY, "test", "tools", "performance.jl"))
    in_timing(f, arguments...; keywords...) =
        Base.invokelatest(getfield(timing, f), arguments...; keywords...)
    @test in_timing(:cpu_list, "0-3,8\n") == [0, 1, 2, 3, 8]
    mktempdir() do directory
        performance, siblings = joinpath(directory, "cpus"), joinpath(directory, "siblings")
        write(performance, "0-15\n")
        write(siblings, "0-1\n")
        @test in_timing(:shared_cpu; taskset = "taskset", performance, siblings, cpus = 32) == 2
        @test in_timing(:shared_cpu; taskset = "taskset", performance = joinpath(directory, "none"),
            siblings = joinpath(directory, "none"), cpus = 4) == 1
        @test in_timing(:shared_cpu; taskset = nothing, performance, siblings, cpus = 32) === nothing
        write(performance, "0-1\n")
        @test in_timing(:shared_cpu; taskset = "taskset", performance, siblings, cpus = 32) === nothing
    end
    @test in_timing(:spread, [100.0, 101.0, 102.0]) ≈ 0.02
    # A stable scenario more than 10 % slower is a slowdown, one 5 to 10 % slower a possible
    # slowdown below the tool's resolution. A round spread above 2 % makes it unstable.
    stable = (0.01, 0.015)
    @test in_timing(:verdict, 0.11, stable) === :slower
    @test in_timing(:verdict, 0.07, stable) === :possible
    @test in_timing(:verdict, 0.05, stable) === :same
    @test in_timing(:verdict, -0.2, stable) === :same
    @test in_timing(:verdict, 0.2, (0.021, 0.01)) === :unstable
    # Exit 1 for a slowdown or a failure, otherwise 2 for a possible slowdown or an
    # unstable scenario, otherwise 0.
    @test in_timing(:status, [:same, :slower, :unstable]) == 1
    @test in_timing(:status, [:same, :failed]) == 1
    @test in_timing(:status, [:same, :possible]) == 2
    @test in_timing(:status, [:same, :unstable]) == 2
    @test in_timing(:status, [:same, :same]) == 0

    # A revision without its own copy uses the working tree's.
    mktempdir() do revision
        corpus = joinpath("test", "support", "scenarios.jl")
        working = joinpath(P.REPOSITORY, corpus)
        marker = "function preservation_corpus("
        @test in_equivalence(:own_copy, revision, corpus; marker) == (; path = working, own = false)
        mkpath(dirname(joinpath(revision, corpus)))
        write(joinpath(revision, corpus), "module CurrentScenarios end\n")
        @test in_equivalence(:own_copy, revision, corpus; marker) == (; path = working, own = false)
        write(joinpath(revision, corpus), "function preservation_corpus(n) end\n")
        @test in_equivalence(:own_copy, revision, corpus; marker) ==
            (; path = joinpath(revision, corpus), own = true)
    end
end
