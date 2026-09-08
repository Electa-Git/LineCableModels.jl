@testitem "Gauntlet / explicit saved benchmarks own comparison direction" tags=[:gauntlet_toolkit] begin
    using JLD2, SHA, TOML
    using LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "docs", "gauntlet_report.jl"))
    using .GauntletArtifacts
    @test occursin("No recorded benchmarks selected", render_gauntlet_report(nothing))
    mktempdir() do root
        frequency = [1.0, 100.0]
        z = reshape(ComplexF64[1, 100], 1, 1, :)
        files = String[]
        function persist(path, document)
            JLD2.save(path, document)
            write(path * ".sha256", bytes2hex(open(sha256, path)) * "  " * basename(path))
        end
        for (name, backend, values) in
            (("a", "fem", z), ("b", "coaxial", reshape(ComplexF64[2, 100], 1, 1, :)),
            ("c", "coaxial", 2z))
            path = joinpath(root, name * ".jld2")
            persist(path,
                Dict{String, Any}("schema_version"=>1, "kind"=>:gauntlet_calculation,
                    "status"=>:complete, "backend"=>backend, "case_id"=>"case_with_underscores",
                    "selection"=>Dict("id"=>name, "earth_impedance"=>name),
                    "formulation"=>(
                        definitions = (earth_impedance = Symbol(name),), options = (;)),
                    "problem"=>Dict("private_input_dump"=>repeat("not for publication", 1000)),
                    "Z"=>values, "Y"=>im .* values, "frequencies"=>frequency,
                    "port_order"=>["core"], "basis"=>:pul, "domain"=>:PhaseDomain,
                    "elapsed_at_completion_seconds"=>99.0, "batch_selection_count"=>2))
            push!(files, path)
        end
        binding(index,
            owner) = Dict("path"=>files[index], "sha256"=>bytes2hex(open(sha256, files[index])), "owner"=>owner)
        entries = [
            Dict("id"=>"forward", "case"=>"case_with_underscores",
                "description"=>"Case_with_underscores <example>",
                "reference"=>binding(1, "external"), "candidate"=>binding(2, "engine")),
            Dict("id"=>"reverse", "case"=>"case_with_underscores",
                "reference"=>binding(2, "engine"), "candidate"=>binding(1, "external")),
            Dict("id"=>"same_backend", "case"=>"case_with_underscores",
                "reference"=>binding(2, "engine"), "candidate"=>binding(3, "engine"))]
        plan = Dict("schema_version"=>1, "collection"=>"fixture", "benchmarks"=>entries,
            "comparison"=>Dict("quantities"=>["Z", "Y", "G"], "bands"=>["all", "dc", "wide"],
                "normalizations"=>["reference_rms", "pointwise"]))
        source = joinpath(root, "benchmarks.toml")
        open(io->TOML.print(io, plan), source, "w")
        original = read.(files)
        output = joinpath(root, "comparisons")
        paths = compare_saved(source; directory = output)
        @test length(paths)==3
        records = read_benchmark.(paths)
        @test records[1]["calculations"].reference.backend=="fem"
        @test records[1]["calculations"].candidate.backend=="coaxial"
        @test records[2]["calculations"].reference.backend=="coaxial"
        @test records[2]["calculations"].candidate.backend=="fem"
        @test records[3]["calculations"].reference.backend==records[3]["calculations"].candidate.backend
        errors = records[1]["reference_comparison"]
        @test only(errors[1].relative) ≈ sqrt(1/10001)
        @test only(errors[2].relative) ≈ sqrt(1/2)
        @test only(records[2]["reference_comparison"][1].relative) ≈ sqrt(1/10004)
        @test records[1]["timings"].reference.scope===:batch_elapsed_at_completion
        @test records[1]["numerical_reference_approval"]===:unreviewed
        before = [(directory, copy(names)) for (directory, _, names) in walkdir(root)]
        summary = render_gauntlet_report(output)
        @test occursin("3 explicitly configured benchmarks", summary)
        @test occursin("Z NRMSE", summary) && occursin("Y NRMSE", summary)
        @test occursin("Z pointwise", summary) && occursin("Y pointwise", summary)
        @test occursin("(1, 1)", summary) && occursin("1=core", summary)
        @test occursin("Case_with_underscores &lt;example&gt;", summary)
        @test occursin("not recorded", summary)
        @test occursin("No stored samples", summary)
        @test length(findall("No stored samples", summary)) == 1
        full_band, slices = split(summary, "### Frequency slices — Z/Y")
        @test occursin("Full-band Z/Y summary", full_band)
        @test !occursin("Band `dc`", full_band)
        @test occursin("Band `dc`", slices) && occursin("Band `wide`", slices)
        @test !occursin("<th>quantity</th>", summary)
        @test !occursin("<th>band</th>", summary)
        @test !occursin("G NRMSE", summary) && !occursin("G pointwise", summary)
        @test any(row -> row.quantity === :G, errors) # Retained, not promoted into the summary.
        selections = Dict(record["calculations"].reference.selection=>index
            for (index, record) in enumerate(records))
        for record in records, operand in record["calculations"]
            get!(selections, operand.selection, length(selections) + 1)
        end
        table = gauntlet_table(records, "all", selections)
        @test nrow(table) == length(records)
        @test names(table) == ["reference", "candidate", "Z NRMSE", "Z pointwise", "Y NRMSE", "Y pointwise"]
        @test startswith(table.reference[1], "fem [")
        @test startswith(table.candidate[1], "coaxial [")
        @test startswith(table.reference[2], "coaxial [")
        @test startswith(table.candidate[2], "fem [")
        @test startswith(table.reference[3], "coaxial [") && startswith(table.candidate[3], "coaxial [")
        @test table[1, "Z NRMSE"] == string(round(100sqrt(1/10001); sigdigits=4), " (1, 1)")
        @test table[1, "Z pointwise"] == string(round(100sqrt(1/2); sigdigits=4), " (1, 1)")
        @test gauntlet_table(records, "wide", selections)[1, "Y NRMSE"] == "missing (no samples)"
        incomplete = deepcopy(first(records))
        filter!(row -> row.quantity === :Z && row.details.normalization === :reference_rms,
            incomplete["reference_comparison"])
        partial = gauntlet_table([incomplete], "all", selections)
        @test partial[1, "Z pointwise"] == "missing (not recorded)"
        @test partial[1, "Y NRMSE"] == "missing (not recorded)"
        @test ncodeunits(summary)<30_000
        for token in
            ("private_input_dump", "not for publication", "<svg", "<img", "case_1.html")
            @test !occursin(token, summary)
        end
        @test read.(files)==original
        @test [(directory, copy(names)) for (directory, _, names) in walkdir(root)]==before
        @test render_gauntlet_report(join((output, output), Sys.iswindows() ? ';' :
                                                            ':'))==summary
        @test_throws r"already exists" compare_saved(source; directory = output)
        delete!(entries[1], "reference")
        open(io->TOML.print(io, plan), source, "w")
        @test_throws r"explicit reference" compare_saved(source; directory = joinpath(root, "invalid"))
        @test !ispath(joinpath(root, "invalid"))
        write(files[1] * ".sha256", "invalid")
        @test_throws r"checksum mismatch" read_calculation(files[1])
        rm(files[1])
        @test_throws r"operand missing" render_gauntlet_report(output)
        @test_throws r"references cannot be inferred" render_gauntlet_report(joinpath(root, "empty") |>
                                                                             mkpath)
    end
end

@testitem "Gauntlet / saved UQ comparisons retain distinct means and deviations" tags=[:gauntlet_toolkit] begin
    using JLD2, SHA, TOML
    include(joinpath(pkgdir(LineCableModels), "docs", "gauntlet_report.jl"))
    using .GauntletArtifacts
    mktempdir() do root
        f=[1.0, 100.0]
        paths=String[]
        for (id, factor) in (("lep", 1.0), ("monte_carlo", 2.0))
            values=(;
                (quantity=>(mean = fill(factor, 1, 1, 2), std = fill(0.1factor, 1, 1, 2))
            for quantity in (:R, :L, :C, :G))...)
            moments=(values = values, frequencies = f, basis = :pul,
                domain = :PhaseDomain, port_order = ["core"])
            path=joinpath(root, "$id.jld2")
            JLD2.jldsave(
                path; schema_version = 1, kind = :gauntlet_moments, status = :complete,
                case_id = "uq_case", backend = :coaxial, problem = (physical = "same",),
                selection = Dict("propagation"=>id), formulation = (
                    definitions = (;), options = (;)),
                frequencies = f, basis = :pul, domain = :PhaseDomain, port_order = ["core"], moments)
            write(path*".sha256", bytes2hex(open(sha256, path)))
            push!(paths, path)
        end
        reference, candidate=read_calculation.(paths)
        source=joinpath(root, "definition.toml")
        write(source, "# fixture")
        benchmark=benchmark_definition(:uq_fixture, :uq_case, :uq, source,
            (id = :uq_case, description = "Mean and standard deviation"),
            benchmark_calculation(:lep, :uq, reference, reference.metadata.formulation),
            benchmark_calculation(:mc, :uq, candidate, candidate.metadata.formulation), UQMomentPolicy(), (;))
        result=compare_saved(benchmark; directory = joinpath(root, "output"))
        record=read_benchmark(result)
        @test record["comparison_policy"].kind===:uq_moments
        @test length(record["reference_comparison"])==8
        @test Set(r.statistic for r in record["reference_comparison"])==Set((:mean, :std))
        @test all(r->only(r.relative)≈1.0, record["reference_comparison"])
        summary=render_gauntlet_report(joinpath(root, "output"))
        @test occursin("mean", summary) && occursin("std", summary)
        @test occursin("### UQ means", summary)
        @test occursin("### UQ standard deviations", summary)
        @test !occursin("Full-band Z/Y summary", summary)
        @test !occursin("pointwise", summary)
    end
end

@testitem "Gauntlet / benchmark Gridspace preserves explicit roles" tags=[:gauntlet_toolkit] begin
    include(joinpath(pkgdir(LineCableModels), "test", "gauntlet", "artifacts.jl"))
    using .GauntletArtifacts
    reference=benchmark_calculation(:fixed_reference, :external, nothing, nothing)
    candidates=Gridspace{BenchmarkCalculation}(
        id->benchmark_calculation(id, :engine, nothing, nothing), (Grid((:one, :two)),))
    benchmarks=Gridspace{OwnedBenchmark}(
        candidate->benchmark_definition(candidate.id, :fixture, :manual, @__FILE__,
            (id = :fixture,), reference, candidate, LineParametersPolicy(), (;)), (Grid(collect(candidates)),))
    materialized=collect(benchmarks)
    @test length(materialized)==2
    @test all(b->b.reference===reference, materialized)
    @test getproperty.(getproperty.(materialized, :candidate), :id)==[:one, :two]
    @test all(b->b.reference.owner===:external, materialized)
end

@testitem "Gauntlet / reporting stays documentation-only" tags=[:gauntlet_toolkit] begin
    using TOML
    root = pkgdir(LineCableModels)
    gauntlet = joinpath(root, "test", "gauntlet")
    include(joinpath(gauntlet, "artifacts.jl"))
    @test !isdefined(GauntletArtifacts, :report)
    @test !isdefined(GauntletArtifacts, :load_report)
    @test isdefined(LineCableModels.Engine, :RMSError)
    @test isdefined(LineCableModels.Engine, :LineParametersBenchmark)
    @test isdefined(LineCableModels.Engine, :compare)
    for file in ("detailed_reports.jl", "reporting.jl", "report.jl", "fem_report.jl",
        "fem_catalogue_report.jl", "plot_fem_matrix_comparison.jl", "plot_fem_admittance_comparison.jl")
        @test !isfile(joinpath(gauntlet, file))
    end
    project = TOML.parsefile(joinpath(gauntlet, "Project.toml"))
    @test !haskey(project["deps"], "CairoMakie")
    source = read(joinpath(root, "docs", "gauntlet_report.jl"), String)
    for token in ("CairoMakie", "GauntletSupport", "load_case",
        "compute(", "cp(", "mkpath(", "write(")
        @test !occursin(token, source)
    end
    @test !occursin("Engine.compare(", source)
    @test !occursin("backend == \"coaxial\"", source)
    @test !occursin("baseline=", source)
    definitions=read(joinpath(gauntlet, "definitions.jl"), String)
    @test !occursin("cannot execute an external", definitions)
    cli = read(joinpath(gauntlet, "cli.jl"), String)
    @test !occursin("gauntlet report", cli)
    @test !occursin("--report", cli)
    @test !occursin("detailed_reports", cli)
end

@testitem "Gauntlet / common benchmark runner accepts either external operand" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport
    using LineCableModels.Engine
    struct StoredFixture
        factor::Float64
    end
    function LineCableModels.Engine.compute(::Nothing, selection::StoredFixture)
        values=fill(complex(selection.factor), 1, 1, 2)
        LineParameters(PhaseDomain, values, im .* values, [1.0, 100.0])
    end
    GauntletSupport._owned_formulation_record(selection::StoredFixture)=(factor = selection.factor,)
    model=load_case(:two_insulated_wires; variation = ExactOverrides(frequencies = [50.0]))
    limits=(Z = (absolute = 10.0, relative = 10.0), Y = (absolute = 10.0, relative = 10.0))
    withenv("LINECABLEMODELS_GAUNTLET_MODE"=>"live") do
        for owners in ((:external, :engine), (:engine, :external))
            benchmark=benchmark_definition(:fixture, model.id, :manual, @__FILE__, model,
                benchmark_calculation(:reference, owners[1], nothing, StoredFixture(1.0)),
                benchmark_calculation(:candidate, owners[2], nothing, StoredFixture(2.0)),
                LineParametersPolicy(normalizations = (:reference_rms, :pointwise)),
                (reference = limits, regression = limits))
            result=run_benchmark(benchmark)
            @test result.metadata.calculations.reference.owner===owners[1]
            @test result.metadata.calculations.candidate.owner===owners[2]
            @test length(result.configured_comparisons)==4
            @test all(row->only(row.error.relative)==1.0, result.configured_comparisons)
        end
    end
end
