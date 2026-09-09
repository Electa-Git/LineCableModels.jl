@testitem "Gauntlet / REPL reports preserve all terms and saved operand identities" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using JLD2, SHA, TOML, DataFrames
    using LineCableModels
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using .GauntletSupport.Gauntlet
    mktempdir() do root
        f = [0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]
        z = cat(([1.0 2.0; 3.0 5.0] .* (1 + im * frequency) for frequency in f)...; dims=3)
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
        files_before = [(dir, copy(names)) for (dir, _, names) in walkdir(root)]
        tables = report(BenchmarkTableDefinition(false), loaded).table
        @test propertynames(tables) == (:calculations, :comparisons, :terms)
        @test nrow(tables.calculations) == 2
        @test tables.calculations.formulation[1].equation === :reference
        @test tables.calculations.formulation[2].equation === :candidate
        @test tables.calculations.selection == [(id=:reference,), (id=:candidate,)]
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
        @test all(reason -> occursin("numerically zero", reason), g.reason)
        @test all(==("S/m"), g.absolute_unit)
        r = filter(row -> row.quantity === :R && row.band === :all, tables.terms)
        @test all(value -> value ≈ 100, r.relative_rms_percent)
        @test all(==("Ω/m"), r.absolute_unit)
        empty_band = filter(row -> row.band == (1e8, 1e9), tables.terms)
        @test all(ismissing, empty_band.absolute_rms)
        @test all(ismissing, empty_band.relative_rms_percent)
        @test all(iszero, empty_band.samples)
        @test all(==(:no_samples), empty_band.status)
        @test tables.comparisons.absolute_rms[1] !== only(loaded.analyses)["reference_comparison"][1].absolute
        @test [(dir, copy(names)) for (dir, _, names) in walkdir(root)] == files_before
        @test all(read(path) == bytes for (path, bytes) in original)
        @test_throws r"explicitly saved" LineCableModels.plot(loaded, (R,); pair=(2, 1))
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
