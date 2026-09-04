@testitem "Quality / observations / one publication and wide tables" tags = [:quality] begin
    root=pkgdir(LineCableModels)
    grammar=read(joinpath(root, "src", "grammar", "observables.jl"), String)
    reports=join(
        read(joinpath(root, "src", "reportbuilder", file), String)
    for file in ("grammar.jl", "tables.jl", "montecarlo.jl", "xlsx.jl")
    )
    makie=join(
        read(path, String)
    for (directory, _, files) in walkdir(
        joinpath(root, "ext", "LineCableModelsMakieExt")
    )
    for file in files if endswith(file, ".jl")
    for path in (joinpath(directory, file),)
    )

    @test parentmodule(LineCableModels.observables) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.ReportBuilder.observation_columns) ===
          LineCableModels.ReportBuilder
    @test occursin("requests::Tuple", grammar)
    @test !occursin("requests::NamedTuple", grammar)
    @test occursin("clip::Bool = true", grammar)

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
    @test occursin("function LineCableModels.plotwindow", makie)
    @test occursin("function _addon_axis!", makie)
    @test occursin(
        "import LineCableModels.PlotBuilder: plot, preview, show_material_scale",
        makie
    )
end
