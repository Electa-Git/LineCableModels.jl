@testitem "Gauntlet / documentation summarizes stored defaults only" tags=[:gauntlet_toolkit] begin
    using JLD2, SHA, TOML
    using LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "docs", "gauntlet_report.jl"))
    @test occursin("No recorded results selected", render_gauntlet_report(nothing))
    mktempdir() do root
        frequency = [0.1, 50.0, 1e6]
        z = reshape(ComplexF64[1 + im, 2 + 2im, 3 + 3im], 1, 1, :)
        y = reshape(ComplexF64[0, im, 2im], 1, 1, :)
        selection = Dict("id"=>"default", "earth_impedance"=>"default")
        variant = Dict("id"=>"nondefault", "earth_impedance"=>"UnknownFormula")
        definitions = Formulation().definitions
        jobs = Dict{String, Any}[]
        documents = Dict{String, Any}[]
        paths = String[]
        function persist(path, document)
            JLD2.save(path, document)
            write(path * ".sha256", bytes2hex(open(sha256, path)) * "  " * basename(path))
        end
        for (backend, factor) in (("coaxial", 1.0), ("fem", 1.1), ("pscad", 0.9))
            job = "case_with_underscores_$backend"
            mkpath(joinpath(root, job))
            push!(jobs, Dict("id"=>job, "case"=>"case_with_underscores", "backend"=>backend,
                "description"=>"Case_with_underscores <example>", "selections"=>[selection, variant]))
            path = joinpath(root, job, "0001.jld2")
            document = Dict{String, Any}("schema_version"=>1, "kind"=>:gauntlet_calculation,
                "status"=>:complete, "backend"=>Symbol(backend), "case_id"=>"case_with_underscores",
                "selection"=>selection, "formulation"=>(; definitions),
                "problem"=>Dict("private_input_dump"=>repeat("not for publication", 1000)),
                "Z"=>factor .* z, "Y"=>factor .* y, "frequencies"=>frequency,
                "port_order"=>["core"], "basis"=>:pul, "domain"=>:PhaseDomain)
            persist(path, document)
            # A nondefault artifact must not even be loaded by the summary.
            write(joinpath(root, job, "0002.jld2"), "not a JLD2 file")
            push!(documents, document)
            push!(paths, path)
        end
        push!(jobs, Dict("id"=>"pending_coaxial", "case"=>"pending", "backend"=>"coaxial",
            "selections"=>[selection]))
        push!(jobs, Dict("id"=>"uq_coaxial", "case"=>"uq", "backend"=>"coaxial",
            "propagation"=>"monte_carlo", "selections"=>[selection]))
        open(joinpath(root, "campaign.toml"), "w") do io
            TOML.print(io, Dict("schema_version"=>1, "jobs"=>jobs))
        end
        before = read.(paths)
        before_files = [(directory, copy(files)) for (directory, _, files) in walkdir(root)]
        defaults = gauntlet_defaults(root)
        @test length(defaults.records) == 3
        @test defaults.pending == [(case="pending", backend="coaxial")]
        comparisons = gauntlet_comparisons(defaults.records)
        @test length(comparisons) == 2
        @test all(row -> row.baseline == "coaxial", comparisons)
        @test Set(row.backend for row in comparisons) == Set(("fem", "pscad"))
        @test isempty(gauntlet_comparisons(filter(record -> record.backend != "coaxial", defaults.records)))
        row = only(filter(row -> row.backend == "fem", comparisons))
        @test row.max_z_percent ≈ 10
        @test row.max_y_percent ≈ 10
        @test compare(defaults.records[1].parameters, defaults.records[2].parameters) isa LineParametersBenchmark
        summary = render_gauntlet_report(root)
        @test occursin("3 completed defaults across 1 cases; 1 default selections unfinished", summary)
        @test occursin("max εZ", summary)
        @test occursin("max εY", summary)
        @test occursin("Case_with_underscores", summary)
        @test occursin("&lt;example&gt;", summary)
        @test ncodeunits(summary) < 12_000
        for token in ("private_input_dump", "not for publication", "<svg", "<img", "<details", "case_1.html")
            @test !occursin(token, summary)
        end
        @test read.(paths) == before
        @test [(directory, copy(files)) for (directory, _, files) in walkdir(root)] == before_files
        @test render_gauntlet_report(join((root, root), Sys.iswindows() ? ';' : ':')) == summary
        original = documents[2]["problem"]
        documents[2]["problem"] = Dict("different_inputs"=>true)
        persist(paths[2], documents[2])
        @test_throws r"inputs or coordinates differ" render_gauntlet_report(root)
        documents[2]["problem"] = original
        persist(paths[2], documents[2])
        documents[2]["formulation"] = (; definitions=merge(definitions, (semicon_admittance=:Ametani2004,)))
        persist(paths[2], documents[2])
        @test length(gauntlet_defaults(root).records) == 2
        documents[2]["formulation"] = (; definitions)
        persist(paths[2], documents[2])
        write(paths[1] * ".sha256", "invalid checksum")
        @test_throws r"checksum mismatch" render_gauntlet_report(root)
        persist(paths[1], documents[1])
        rm(paths[1] * ".sha256")
        @test length(gauntlet_defaults(root).pending) == 2
    end
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
    for token in ("CairoMakie", "GauntletSupport", "load_case", "compute(", "cp(", "mkpath(", "write(")
        @test !occursin(token, source)
    end
    cli = read(joinpath(gauntlet, "cli.jl"), String)
    @test !occursin("gauntlet report", cli)
    @test !occursin("--report", cli)
    @test !occursin("detailed_reports", cli)
end
