@testitem "Core / architecture / removed infrastructure stays removed" tags=[:unit] begin
    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    source_files=filter(endswith(".jl"),
        [joinpath(directory, file)
         for (directory, _, files) in walkdir(source_root)
         for file in files])
    source=Dict(relpath(path, root)=>read(path, String) for path in source_files)

    @test !isdefined(LineCableModels, :Commons)
    @test !isdefined(LineCableModels, :Utils)
    @test !ispath(joinpath(source_root, "commons"))
    @test !ispath(joinpath(source_root, "utils"))
    @test isempty(filter(path -> basename(path) == "helpers.jl", source_files))

    forbidden=(
        "resolve_T", "coerce_to_T", "BASE_FLOAT", "REALSCALAR", "COMPLEXSCALAR",
        "@construct", "@parameterize", "@measurify", "Meta.parse", "global_logger(",
        "line_component_quantity", "line_component_unit", "parse_kwargs",
        "resolve_input", "axis_quantity", "axis_unit", "series_data",
        "ComponentMetadata", "UnitSpec", "DEFAULT_QUANTITY_UNITS",
        "_mc_selector", "AbstractObservable", "ObservableDescriptor",
        "ObservationRequest", "AbstractDetails", "DetailsRegistry"
    )
    for token in forbidden
        @test all(!occursin(token, contents) for contents in values(source))
    end
    @test all(!occursin("@assert", contents) for contents in values(source))

    macro_path=joinpath("src", "parametricbuilder", "macros.jl")
    macro_facades=Set((
        joinpath("src", "LineCableModels.jl"),
        joinpath("src", "parametricbuilder", "ParametricBuilder.jl")
    ))
    for (path, contents) in source
        path == macro_path && continue
        if path in macro_facades
            @test occursin("@relax", contents)
            @test !occursin(r"\brecast\b", contents)
            continue
        end
        @test !occursin("@relax", contents)
        @test !occursin(r"\brecast\b", contents)
    end

    maintained=filter(pair -> pair.first != "src/retired.jl", collect(source))
    @test all(!occursin("calc_", contents) for (_, contents) in maintained)

    eager_files=(
        joinpath("src", "datamodel", "circstrands.jl"),
        joinpath("src", "datamodel", "rectstrands.jl"),
        joinpath("src", "datamodel", "strip.jl"),
        joinpath("src", "datamodel", "tubular.jl"),
        joinpath("src", "datamodel", "insulator.jl"),
        joinpath("src", "datamodel", "semicon.jl"),
        joinpath("src", "datamodel", "conductorgroup.jl"),
        joinpath("src", "datamodel", "insulatorgroup.jl"),
        joinpath("src", "datamodel", "cablecomponent", "cablecomponent.jl"),
        joinpath("src", "datamodel", "cabledesign", "cabledesign.jl")
    )
    for path in eager_files
        contents=source[path]
        @test !occursin(r"\btemperature\s*::", contents)
        @test !occursin(r";\s*temperature\s*=", contents)
    end
end

@testitem "Core / architecture / ownership-centered recursive layout" tags=[:unit] begin
    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    extension_root=joinpath(root, "ext")
    source_files=filter(endswith(".jl"),
        [joinpath(directory, file)
         for (directory, _, files) in walkdir(source_root)
         for file in files])
    source=Dict(relpath(path, root)=>read(path, String) for path in source_files)

    expected_entries=(
        joinpath("src", "grammar", "Grammar.jl"),
        joinpath("src", "units", "Units.jl"),
        joinpath("src", "validation", "Validation.jl"),
        joinpath("src", "plotbuilder", "PlotBuilder.jl"),
        joinpath("src", "datamodel", "DataModel.jl"),
        joinpath("src", "engine", "Engine.jl"),
        joinpath("src", "uq", "UQ.jl"),
        joinpath("src", "reportbuilder", "ReportBuilder.jl"),
        joinpath("src", "importexport", "ImportExport.jl"),
        joinpath("ext", "LineCableModelsMakieExt", "LineCableModelsMakieExt.jl")
    )
    @test all(isfile(joinpath(root, path)) for path in expected_entries)

    expected_owned_files=(
        joinpath("src", "datamodel", "packing.jl"),
        joinpath("src", "datamodel", "cabledesign", "cabledesign.jl"),
        joinpath("src", "datamodel", "preview", "cable.jl"),
        joinpath("src", "engine", "lineparameters", "lineparameters.jl"),
        joinpath("src", "engine", "lineparameters", "quantities.jl"),
        joinpath("src", "engine", "earthreturn.jl"),
        joinpath("src", "uq", "montecarlo", "compute.jl"),
        joinpath("src", "reportbuilder", "grammar.jl"),
        joinpath("src", "reportbuilder", "tables.jl"),
        joinpath("src", "reportbuilder", "xlsx.jl"),
        joinpath("src", "importexport", "pscad", "pscad.jl"),
        joinpath("ext", "LineCableModelsMakieExt", "UIComponents.jl"),
        joinpath("dev", "plotting", "Project.toml")
    )
    @test all(isfile(joinpath(root, path)) for path in expected_owned_files)

    obsolete_paths=(
        joinpath("src", "Grammar.jl"),
        joinpath("src", "unithandler"),
        joinpath("src", "plotbuilder", "backendhandler"),
        joinpath("src", "plotbuilder", "uicomponents"),
        joinpath("src", "datamodel", "strands_handler.jl"),
        joinpath("src", "engine", "solver.jl"),
        joinpath("src", "uq", "plotspecs.jl"),
        joinpath("src", "engine", "lineparameters", "dataframe.jl"),
        joinpath("src", "uq", "dataframe.jl"),
        joinpath("src", "importexport", "xlsx.jl"),
        joinpath("integration", "plotting")
    )
    @test all(!ispath(joinpath(root, path)) for path in obsolete_paths)

    @test !isdefined(LineCableModels, :UnitHandler)
    @test !isdefined(LineCableModels.PlotBuilder, :BackendHandler)
    @test isdefined(LineCableModels.DataModel, :PhysicalFillLimit)
    @test !isdefined(LineCableModels.Validation, :PhysicalFillLimit)

    @test parentmodule(LineCableModels.description) === LineCableModels.Engine
    @test parentmodule(LineCableModels.Z) === LineCableModels.Engine
    @test parentmodule(LineCableModels.frequencies) === LineCableModels.Engine
    @test parentmodule(LineCableModels.ncables) === LineCableModels.DataModel
    @test parentmodule(LineCableModels.nphases) === LineCableModels.DataModel
    @test parentmodule(LineCableModels.compute) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.observe) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.observables) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.Units.quantity) === LineCableModels.Units
    @test parentmodule(LineCableModels.report) === LineCableModels.ReportBuilder
    @test parentmodule(LineCableModels.ReportArtifact) === LineCableModels.ReportBuilder
    @test parentmodule(LineCableModels.XLSXReport) === LineCableModels.ReportBuilder
    @test parentmodule(LineCableModels.validate) === LineCableModels.Validation

    report_grammar=source[joinpath("src", "reportbuilder", "grammar.jl")]
    @test !occursin(r"(?:encode|write)\(::AbstractReportDefinition", report_grammar)
    importexport_index=source[joinpath("src", "importexport", "ImportExport.jl")]
    reportbuilder_index=source[joinpath("src", "reportbuilder", "ReportBuilder.jl")]
    @test !occursin("using XLSX", importexport_index)
    @test !occursin("include(\"xlsx.jl\")", importexport_index)
    @test occursin("using XLSX", reportbuilder_index)
    @test occursin("include(\"xlsx.jl\")", reportbuilder_index)

    for pattern in (
        r"details\s*::\s*Dict",
        r"details\s*::\s*Any",
        r"function\s+computation_details\s*\([^,)]*\)"
    )
        @test all(!occursin(pattern, contents) for contents in values(source))
    end

    presentation_paths=filter(keys(source)) do path
        startswith(path, joinpath("src", "plotbuilder")) ||
            startswith(path, joinpath("src", "reportbuilder")) ||
            endswith(path, "plotdefinition.jl") ||
            endswith(path, "comparisonplot.jl") ||
            endswith(path, joinpath("montecarlo", "plot.jl"))
    end
    for path in presentation_paths
        contents=source[path]
        @test !occursin(r"\.(Z|Y|f)\b", contents)
        @test !occursin(r"\b(Ω/m|Ω/km|H/m|mH/km|S/m|S/km|F/m|μF/km)\b", contents)
    end

    @test all(!occursin("@eval", contents) for contents in values(source))
    @test all(!occursin(r"^\s*(using|import)\s+.*Makie"m, contents)
    for contents in values(source))
    @test all(!occursin("@eval", read(path, String))
    for (directory, _, files) in walkdir(extension_root)
    for file in files if endswith(file, ".jl")
    for path in (joinpath(directory, file),))
end

@testitem "ParametricBuilder / relax / isolated future provision" tags=[:unit] begin
    using LineCableModels.ParametricBuilder: @relax

    @relax struct FutureValue{T <: Real}
        scalar::T
        values::Vector{T}
        label::Symbol
    end

    value=FutureValue(1.0f0, [2.0, 3.0], :future)
    @test value isa FutureValue{Float64}
    @test value.scalar == 1.0
    @test value.values == [2.0, 3.0]
    converted=convert(FutureValue{Float32}, value)
    @test converted isa FutureValue{Float32}
    @test converted.values == Float32[2, 3]
    @test converted.label === :future
end
