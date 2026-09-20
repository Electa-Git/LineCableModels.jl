
@testitem "Gauntlet / explicit saved benchmarks own comparison direction" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using JLD2, SHA, TOML, Measurements
    using LineCableModels.Engine
    using .GauntletSupport: Gauntlet
    include(joinpath(pkgdir(LineCableModels), "docs", "gauntlet_report.jl"))
    using .GauntletSupport.Gauntlet
    @test occursin("No published benchmark artifacts", render_gauntlet_report(nothing))
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
        @test !haskey(records[1],"numerical_reference_approval")
        before = [(directory, copy(names)) for (directory, _, names) in walkdir(root)]
        summary = render_gauntlet_report(output)
        @test occursin("3 complete benchmarks",summary)
        @test all(band -> occursin(band,summary),("all","dc","wide"))
        @test !occursin("<svg",summary) && !occursin("data:image",summary)
        @test !occursin("<img",summary) && !occursin("![",summary)
        @test !occursin("Absolute RMS [",summary)
        @test !occursin("Worst terms and comparison counts",summary)
        @test !occursin("Scientific formula descriptions",summary)
        @test !isdefined(@__MODULE__,:gauntlet_table)
        tables=[report(LineCableModels.ReportBuilder.BenchmarkTableDefinition(),read_benchmark(path;load_results=true)).tables for path in paths]
        row=only(filter(row -> row.quantity===:Z && row.band===:all && row.normalization===:reference_rms,tables[1].maxima))
        @test row.maximum_relative_rms_percent ≈ 100sqrt(1/10001)
        reverse_row=only(filter(row -> row.quantity===:Z && row.band===:all && row.normalization===:reference_rms,tables[2].maxima))
        @test reverse_row.maximum_relative_rms_percent ≈ 100sqrt(1/10004)
        empty_rows=filter(row -> row.band===:wide,tables[1].maxima)
        @test all(ismissing,empty_rows.maximum_relative_rms_percent)
        @test all(iszero,empty_rows.samples)
        @test all(row -> row.unavailable == row.term_count,eachrow(empty_rows))
        @test all(table -> Set(table.terms.quantity)==Set((:Z,:Y,:G)),tables)
        # Publication must embed the owner's complete numeric table, not the
        # old quantities-as-columns text summary or a second RMS aggregation.
        for table in tables,feature in table.features
            rendered=sprint((io,frame) -> show(IOContext(io,:limit=>false),MIME"text/html"(),
                frame;summary=false,eltypes=false),feature.relative)
            @test occursin(rendered,summary)
        end
        for token in
            ("private_input_dump", "not for publication", "all slots :default", "remaining slots :default")
            @test !occursin(token, summary)
        end
        @test read.(files)==original
        @test [(directory, copy(names)) for (directory, _, names) in walkdir(root)]==before
        @test render_gauntlet_report(join((output, output), Sys.iswindows() ? ';' :
                                                            ':'))==summary
        @test all(a != b for (a,b) in zip(compare_saved(source; directory = output),paths))
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


@testitem "Gauntlet / saved UQ comparisons retain distinct means and deviations" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using JLD2, SHA, TOML, Measurements
    using .GauntletSupport: Gauntlet
    include(joinpath(pkgdir(LineCableModels), "docs", "gauntlet_report.jl"))
    using .GauntletSupport.Gauntlet
    mktempdir() do root
        f=[1.0, 100.0]
        paths=String[]
        for (id,factor) in (("reference",1.),("candidate",2.))
            value=measurement(factor,.1factor^2)
            omega=reshape(2pi.*f,1,1,:)
            core=LineParameters(fill(value,1,1,2).+im.*omega.*value,
                fill(value,1,1,2).+im.*omega.*value,f)
            result=LinearErrorResult(LinearError(Formulation()),[core])
            path=joinpath(root,"$id.jld2")
            JLD2.jldsave(path;schema_version=3,kind=:gauntlet_uncertainty,status=:complete,
                case_id="uq_case",backend=:coaxial,problem=(physical="same",),
                selection=(propagation=:linear_error,),formulation=NamedTuple(result).formulation,
                frequencies=f,basis=:pul,domain=:PhaseDomain,port_order=["core"],
                scientific_result=LineCableModels.ImportExport.serialize_value(result))
            write(path*".sha256",bytes2hex(open(sha256,path)))
            push!(paths,path)
        end
        reference,candidate=read_calculation.(paths)
        @test reference.metadata.recovery===:portable
        @test candidate.metadata.recovery===:portable
        @test reference.result isa LinearErrorResult
        source=joinpath(root, "definition.toml")
        write(source, "# fixture")
        benchmark=benchmark_definition(:uq_fixture, :uq_case, :uq, source,
            (id = :uq_case, description = "Mean and standard deviation"),
            BenchmarkCalculation(:reference, reference, reference.metadata.formulation),
            BenchmarkCalculation(:candidate, candidate, candidate.metadata.formulation), (quantities=(:R, :L, :C, :G), statistics=(:mean, :std)), (;))
        result=compare_saved(benchmark; directory = joinpath(root, "output"))
        record=read_benchmark(result)
        @test length(record["comparison_settings"].requests)==8
        @test length(record["reference_comparison"])==40
        @test Set(r.statistic for r in record["reference_comparison"])==Set((:mean, :std))
        # Mean and spread deliberately disagree by different amounts: swapping
        # them must fail even when both have the same matrix and band coordinates.
        @test all(r->only(r.relative)≈(r.statistic===:mean ? 1.0 : 3.0),
            filter(row -> row.details.band===:all,record["reference_comparison"]))
        loaded=read_benchmark(result;load_results=true)
        tables=report(LineCableModels.ReportBuilder.BenchmarkTableDefinition(false),loaded).tables
        @test Set(tables.comparisons.statistic)==Set((:mean,:std))
        @test length(tables.comparisons.quantity)==40
        @test all(row -> only(row.relative_rms_percent)≈(row.statistic===:mean ? 100 : 300),
            eachrow(filter(row -> row.band===:all,tables.comparisons)))
        for feature in tables.features
            @test only(feature.relative.all)≈(feature.statistic===:mean ? 100 : 300)
        end
        summary=render_gauntlet_report(joinpath(root, "output"))
        @test occursin("mean", summary) && occursin("std", summary)
        @test occursin("Relative RMS", summary)
        @test !occursin("<svg",summary)
        @test !occursin("Full-band comparisons", summary)
        @test !occursin("pointwise", summary)
        @test !occursin("Retained UQ statistics",summary)
        @test !occursin("mean_standard_error",summary)
        @test !occursin("configuration 1",summary)
        # First-order results without sampling evidence must not acquire an MC
        # trial count or a fabricated CDF bound simply because they are reported.
        @test !occursin("MC sampling workload",summary)
        @test !occursin("CDF precision",summary)
    end
end


@testitem "Gauntlet / benchmark Gridspace preserves explicit roles" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    reference=BenchmarkCalculation(:fixed_reference, nothing, nothing)
    candidates=Gridspace{BenchmarkCalculation}(
        id->BenchmarkCalculation(id, nothing, nothing), (Grid((:one, :two)),))
    benchmarks=Gridspace{BenchmarkDefinition}(
        candidate->benchmark_definition(candidate.id, :fixture, :manual, @__FILE__,
            (id = :fixture,), reference, candidate, (;), (;)), (Grid(collect(candidates)),))
    materialized=collect(benchmarks)
    @test length(materialized)==2
    @test all(b->b.reference===reference, materialized)
    @test getproperty.(getproperty.(materialized, :candidate), :id)==[:one, :two]
    @test all(b->!hasproperty(b.reference,:owner), materialized)
end
