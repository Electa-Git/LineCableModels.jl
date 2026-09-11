@testitem "Quality / observations / one publication and wide tables" tags = [:quality] begin
    root=pkgdir(LineCableModels)
    reports=join(
        read(joinpath(root, "src", "reportbuilder", file), String)
    for file in ("grammar.jl", "tables.jl", "montecarlo.jl", "xlsx.jl")
    )
    @test parentmodule(LineCableModels.observables) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.ReportBuilder.observation_columns) ===
          LineCableModels.ReportBuilder
    # The dispatch guard in orchestrators.jl and the observation behavior tests
    # cover request types and clipping; source spelling is not a second API.

    for token in (
        "zero_atol",
        "relative_status",
        "rms_absolute",
        "rms_relative",
        "line_requests",
        "line_parent"
    )
        @test !occursin(token, reports)
    end
    @test !occursin(r"\bsource\.(?:stats|sample_values|histogram_values)\b", reports)
    @test !occursin(r"\b(?:family|quantity|value|unit)\s*=\s*[^=]",
        read(
            joinpath(root, "src", "reportbuilder", "tables.jl"), String
        ))
    @test occursin("observation_columns(table)", reports)

    @test !occursin("Figure(", read(
        joinpath(root, "ext", "LineCableModelsMakieExt", "montecarlo.jl"), String
    ))
    # Native window/addon behavior is exercised in the visual environment.
    # Private helper names and the spelling/order of imports are not contracts.
end
