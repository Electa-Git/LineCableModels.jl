@testitem "Gauntlet / REPL reports preserve all terms and saved operand identities" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using JLD2, SHA, TOML, DataFrames
    using LineCableModels
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using .GauntletSupport.Gauntlet
    mktempdir() do root
        f = [0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]
        z = [complex(11i+7j+k,2i+5j+3k) for i in 1:2,j in 1:2,k in eachindex(f)]
        y = im .* z .- real.(im .* z) # exactly zero G
        original = Dict{String, Vector{UInt8}}()
        operands = map(("reference", "candidate"), (1.0, 2.0)) do role, factor
            path = joinpath(root, "$role.jld2")
            JLD2.jldsave(path; schema_version=1, kind=:gauntlet_calculation,
                status=:complete, backend="coaxial", case_id="matrix_case",
                selection=(id=Symbol(role),), formulation=(equation=Symbol(role), controls=(rtol=1e-8,)),
                problem=(identity="same",), Z=factor .* z,
                Y=role == "reference" ? y : factor .* y .+ 1e-13,
                frequencies=f, port_order=["first", "second"], basis=:pul, domain=:PhaseDomain)
            digest = bytes2hex(open(sha256, path))
            write(path * ".sha256", digest)
            original[path] = read(path)
            Dict("path"=>path, "sha256"=>digest)
        end
        source = joinpath(root, "definition.toml")
        open(source, "w") do io
            TOML.print(io, Dict("schema_version"=>1, "collection"=>"repl",
                "comparison"=>Dict("quantities"=>["Z", "Y", "R", "L", "G", "C"],
                    "normalizations"=>["reference_rms", "pointwise"],
                    "bands"=>["all", [1e8, 1e9]]),
                "benchmarks"=>[Dict("id"=>"matrix_report", "case"=>"matrix_case",
                    "reference"=>operands[1], "candidate"=>operands[2])]))
        end
        snapshot = only(compare_saved(source; directory=joinpath(root, "analysis")))
        loaded = read_benchmark(snapshot; load_results=true)
        # Reanalysis identity includes scientific semantics, not runtime sources.
        @test only(compare_saved(source; directory=joinpath(root, "analysis"))) == snapshot
        current = report(BenchmarkTableDefinition(; only(loaded.analyses)["comparison_settings"]...),
            (reference=loaded.reference, candidate=loaded.candidate))
        operands_for_definition = map((:reference, :candidate)) do role
            operand = getproperty(loaded, role)
            Gauntlet.BenchmarkCalculation(role, operand, operand.metadata.formulation)
        end
        definition = Gauntlet.benchmark_definition(:matrix_report, :matrix_case, :repl,
            source, (id=:matrix_case, description="matrix_case"), operands_for_definition...,
            current.published.settings, (;))
        @test Gauntlet.record_benchmark(definition, current;
            directory=joinpath(root, "analysis")) == snapshot
        @test all(==(LineCableModels.Engine.OBSERVABLE_RESOLUTION_REVISION), current.table.terms.resolution_revision)
        @test all(read(path) == bytes for (path, bytes) in original)
        files_before = [(dir, copy(names)) for (dir, _, names) in walkdir(root)]
        tables = report(BenchmarkTableDefinition(false), loaded).table
        @test propertynames(tables) == (:calculations,:formulations,:formula_details,:comparisons,:terms,:maxima,:summary,:features,
            :execution,:source_timings,:performance,:performance_samples,:performance_environment,
            :performance_policy,:performance_comparison,:statistics,:sampling,:mean_sampling_precision,:overview)
        @test nrow(tables.calculations) == 2
        @test loaded.reference.metadata.formulation.equation === :reference
        @test loaded.candidate.metadata.formulation.equation === :candidate
        @test loaded.reference.metadata.selection == (id=:reference,)
        @test loaded.candidate.metadata.selection == (id=:candidate,)
        @test all(column -> all(value -> value isa Union{Number,Symbol,AbstractString,Missing},column),
            eachcol(tables.calculations))
        @test nrow(tables.comparisons) == 24
        @test nrow(tables.terms) == 96
        @test Set(tables.terms.quantity) == Set((:Z, :Y, :R, :L, :G, :C))
        @test Set(zip(tables.terms.row, tables.terms.column)) == Set(((1, 1), (1, 2), (2, 1), (2, 2)))
        for (index, saved) in enumerate(only(loaded.analyses)["reference_comparison"])
            @test isequal(tables.comparisons.absolute_rms[index], saved.absolute)
            @test isequal(tables.comparisons.relative_rms_percent[index], 100 .* saved.relative)
            @test tables.comparisons.port_order[index] == ["first", "second"]
        end
        g = filter(row -> row.quantity === :G && row.band === :all, tables.terms)
        @test all(ismissing, g.relative_rms_percent)
        @test all(>(0), g.absolute_rms)
        @test all(==(:reference_below_tolerance), g.status)
        @test all(reason -> occursin("below declared resolution", reason), g.reason)
        @test all(==("S/m"), g.absolute_unit)
        r = filter(row -> row.quantity === :R && row.band === :all, tables.terms)
        @test all(value -> value ≈ 100, r.relative_rms_percent)
        @test all(==("Ω/m"), r.absolute_unit)
        empty_band = filter(row -> row.band == string((1e8, 1e9)), tables.terms)
        @test !isempty(empty_band)
        @test all(ismissing, empty_band.absolute_rms)
        @test all(ismissing, empty_band.relative_rms_percent)
        @test all(iszero, empty_band.samples)
        @test all(==(:no_samples), empty_band.status)
        @test tables.comparisons.absolute_rms[1] !== only(loaded.analyses)["reference_comparison"][1].absolute
        @test [(dir, copy(names)) for (dir, _, names) in walkdir(root)] == files_before
        @test all(read(path) == bytes for (path, bytes) in original)
        default_source=joinpath(root,"default_request.toml")
        default_request=TOML.parsefile(source)
        delete!(default_request,"comparison")
        open(io -> TOML.print(io,default_request),default_source,"w")
        default_snapshot=only(compare_saved(default_source;directory=joinpath(root,"default_analysis")))
        selected=read_benchmark(default_snapshot)
        @test selected["comparison_settings"] == BenchmarkTableDefinition().settings
        @test length(selected["reference_comparison"]) == 30
        @test_throws r"Plotting is optional" LineCableModels.plot(loaded, (R,))
        modified = deepcopy(loaded)
        modified.reference.result.Z.values[1, 2, 1] += 1
        @test_throws r"modified after loading" report(BenchmarkTableDefinition(false), modified)
        changed = deepcopy(loaded)
        entry = only(changed.analyses)
        entry["calculations"] = merge(entry["calculations"],
            (reference=merge(entry["calculations"].reference, (sha256="wrong",)),))
        @test_throws r"differs from the loaded operand" report(BenchmarkTableDefinition(false), changed)
        # Relative bindings remain correct when the entire retained tree moves.
        parent = mktempdir()
        try
            moved = joinpath(parent, "moved")
            cp(root, moved)
            moved_snapshot = joinpath(moved, relpath(snapshot, root))
            relocated = read_benchmark(moved_snapshot; load_results=true)
            @test relocated.reference.metadata.path == joinpath(moved, "reference.jld2")
            @test isequal(report(BenchmarkTableDefinition(false), relocated).table.terms, tables.terms)
        finally
            rm(parent; recursive=true)
        end
    end
end

@testitem "Gauntlet / saved performance is bound to workloads and its original session" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using JLD2,SHA
    using .GauntletSupport.Gauntlet
    mktempdir() do root
        path=joinpath(root,"performance.jld2")
        calculations=(reference=(problem=:a,formulation=:mc),candidate=(problem=:a,formulation=:lep))
        performance=(reference=(calculation=calculations.reference,median_seconds=3.0),
            candidate=(calculation=calculations.candidate,median_seconds=1.0))
        session=(id="original",started_at="2026-09-13T12:00:00")
        JLD2.jldsave(path;schema_version=1,performance,session)
        write(path*".sha256",bytes2hex(open(sha256,path)))
        original=read(path)
        retained=Gauntlet.read_benchmark(path,Val(:performance);calculations)
        @test retained.checksum_verified===true
        @test retained.workload_verified===true
        @test retained.session==session
        @test read(path)==original
        @test_throws r"workload differs" Gauntlet.read_benchmark(path,Val(:performance);
            calculations=merge(calculations,(candidate=(problem=:b,formulation=:lep),)))
        write(path*".sha256","wrong")
        @test_throws r"checksum mismatch" Gauntlet.read_benchmark(path,Val(:performance);calculations)
    end
end
